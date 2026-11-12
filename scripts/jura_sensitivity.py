"""Independently verify scores, summarise selection and assess target sensitivity."""
from pathlib import Path
import json
import os
import subprocess
import sys
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.stats import wasserstein_distance
from sklearn.cluster import KMeans

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT/'build/jura-export'
SUMMARY = ROOT/'revision/results/supplementary/jura'
RESULTS.mkdir(parents=True, exist_ok=True)
SUMMARY.mkdir(parents=True, exist_ok=True)
KEY = ['dataset','target','repeat_id','model']
GROUP = ['dataset','target','model']
SCHEMES = ['random','block','buffered','knndm']

def export_details():
    script = r'''
source('revision/R/local_libraries.R')
root <- Sys.getenv('APPLICATION_RESULT_ROOT')
bundle <- readRDS('revision/results/application-jura.rds')
for (name in c('results','predictions','fits')) write.csv(bundle[[name]],file.path(root,paste0(name,'.csv')),row.names=FALSE)
files <- bundle$details
out <- lapply(files,function(p) {
  x <- p; g <- x$geometry
  list(dataset=x$results$dataset[1],target=x$results$target[1],repeat_id=x$results$repeat_id[1],
       train_rows=g$train_rows,target_rows=g$target_rows,xy=g$xy,target_xy=g$target_xy,
       folds=g$folds,dc=g$dc,dp=g$dp,buffer_km=g$buffer_km,knndm_warnings=g$knndm_warnings,
       variograms=lapply(x$full_models,function(f) list(empirical=f$dependence$variogram,
         range_km=f$dependence$range,psill=f$dependence$psill,nugget=f$dependence$nugget)))
})
jsonlite::write_json(out,file.path(root,'geometry.json'),auto_unbox=TRUE,digits=16)
'''
    subprocess.run(['Rscript','--vanilla','-e',script],check=True,
                   cwd=ROOT,env={**os.environ,'TZ':'UTC','APPLICATION_RESULT_ROOT':str(RESULTS)})

def selection_weights(frame,rule):
    if rule.startswith('fixed_'):
        return (frame['scheme'].to_numpy()==rule[6:]).astype(float)
    if rule=='best_observed':
        values=frame['absolute_error'].to_numpy()
    else:
        values=frame[rule].to_numpy()
        if rule in ('A','J1','J2','J2_residual'):
            values=np.abs(values)
    selected=values<=values.min()+1e-12
    return selected.astype(float)/selected.sum()

def main():
    export_details()
    raw=pd.read_csv(RESULTS/'results.csv')
    pred=pd.read_csv(RESULTS/'predictions.csv')
    fits=pd.read_csv(RESULTS/'fits.csv')
    details=json.loads((RESULTS/'geometry.json').read_text())
    assert len(details)==120 and len(raw)==960 and len(fits)==5040
    assert not raw.duplicated(KEY+['scheme']).any()
    assert np.isfinite(raw.select_dtypes('number')).all().all()
    assert np.allclose(pred.squared_error,(pred.observed-pred.predicted)**2,rtol=1e-12,atol=1e-12)
    cv=pred[pred.role=='cv'].groupby(KEY+['scheme']).squared_error.mean()
    target=pred[pred.role=='target'].groupby(KEY).squared_error.mean()
    maximum_score_difference=0.
    for r in raw.itertuples(index=False):
        key=(r.dataset,r.target,r.repeat_id,r.model)
        errs=[abs(r.cv_mse-cv.loc[key+(r.scheme,)]),abs(r.target_mse-target.loc[key]),
              abs(r.gap-(r.target_mse-r.cv_mse)),abs(r.absolute_error-abs(r.gap))]
        maximum_score_difference=max(maximum_score_difference,*errs)
    assert maximum_score_difference<1e-10
    max_distance_error=max_w_error=max_j2_error=0.
    for g in details:
        xy=np.asarray(g['xy']); q=np.asarray(g['target_xy'])
        dd=cdist(xy,xy); dq=cdist(q,xy).min(axis=1)
        assert set(g['train_rows']).isdisjoint(g['target_rows'])
        max_distance_error=max(max_distance_error,float(np.max(np.abs(dq-g['dp']))))
        for scheme,folds in g['folds'].items():
            dc=np.full(len(xy),np.nan)
            indices=[]
            for fold in folds:
                tr=np.atleast_1d(fold['train']).astype(int)-1
                te=np.atleast_1d(fold['test']).astype(int)-1
                assert not set(tr)&set(te)
                indices.extend(te.tolist())
                dc[te]=dd[np.ix_(te,tr)].min(axis=1)
                if scheme=='buffered':
                    assert dd[np.ix_(te,tr)].min()>g['buffer_km']-1e-8
            assert sorted(indices)==list(range(len(xy))) and np.isfinite(dc).all()
            max_distance_error=max(max_distance_error,float(np.max(np.abs(dc-g['dc'][scheme]))))
            row=raw[(raw.dataset==g['dataset'])&(raw.target==g['target'])&
                    (raw.repeat_id==g['repeat_id'])&(raw.model=='ordinary')&(raw.scheme==scheme)].iloc[0]
            max_w_error=max(max_w_error,abs(wasserstein_distance(dc,dq)-row.W))
            j2=row.signal**2*(np.mean(np.exp(-2*dc/row.range_km))-np.mean(np.exp(-2*dq/row.range_km)))
            max_j2_error=max(max_j2_error,abs(j2-row.J2))
    assert max_distance_error<1e-7 and max_w_error<1e-7 and max_j2_error<1e-9
    # Fixed target models must remain identical while only CV partitions vary.
    fixed=pred[(pred.role=='target')&(pred.target!='uniform')]
    variation=fixed.groupby(['dataset','target','model','site_id']).predicted.agg(lambda x:x.max()-x.min()).max()
    assert variation<1e-12
    checks={'configuration_count':len(details),'candidate_evaluations':len(raw),
            'covariance_fits':len(fits),'max_score_reconstruction_error':maximum_score_difference,
            'max_distance_reconstruction_error_km':max_distance_error,
            'max_W_reconstruction_error_km':max_w_error,'max_J2_reconstruction_error':max_j2_error,
            'max_fixed_target_prediction_change':float(variation),
            'fold_coverage_and_target_separation':'passed',
            'finite_predictions_and_scores':'passed'}
    (SUMMARY/'verification.json').write_text(json.dumps(checks,indent=2)+'\n')
    baselines=[]
    for g in details:
        dataset=g['dataset']
        d=pd.read_csv(ROOT/'revision/data/jura-sites.csv')
        y=d.iloc[np.asarray(g['train_rows'])-1].response.to_numpy()
        yq=d.iloc[np.asarray(g['target_rows'])-1].response.to_numpy()
        baseline=float(np.mean((yq-y.mean())**2))
        frame=raw[(raw.dataset==dataset)&(raw.target==g['target'])&(raw.repeat_id==g['repeat_id'])]
        for model,f in frame.groupby('model'):
            baselines.append({'dataset':dataset,'target':g['target'],'repeat_id':g['repeat_id'],
                'model':model,'training_mean_target_mse':baseline,'model_target_mse':f.target_mse.iloc[0],
                'relative_gain_vs_training_mean':1-f.target_mse.iloc[0]/baseline})
    baseline_frame=pd.DataFrame(baselines)
    baseline_frame.to_csv(SUMMARY/'target-baseline-comparison.csv',index=False)
    baseline_frame.groupby(GROUP).agg(mean_baseline_mse=('training_mean_target_mse','mean'),
        mean_model_mse=('model_target_mse','mean'),mean_relative_gain=('relative_gain_vs_training_mean','mean')).reset_index().to_csv(
        SUMMARY/'target-baseline-summary.csv',index=False)
    rules=['W','J2','A','J1','U1','U2','J2_residual']+['fixed_'+s for s in SCHEMES]+['best_observed']
    selected=[]; weights=[]
    for key,frame in raw.groupby(KEY,sort=True):
        for rule in rules:
            w=selection_weights(frame,rule)
            selected.append(dict(zip(KEY,key))|{'rule':rule,'absolute_error':float(w@frame.absolute_error),
                'cv_mse':float(w@frame.cv_mse),'target_mse':float(frame.target_mse.iloc[0]),
                'gap':float(w@frame.gap),'relative_error':float(w@frame.absolute_error/frame.target_mse.iloc[0]),
                'n_tied':int((w>0).sum()),'selected_schemes':'+'.join(frame.loc[w>0,'scheme'])})
            for scheme,weight in zip(frame.scheme,w):
                weights.append(dict(zip(KEY,key))|{'rule':rule,'scheme':scheme,'weight':float(weight)})
    selected=pd.DataFrame(selected); weights=pd.DataFrame(weights)
    selected.to_csv(SUMMARY/'selection-evaluations.csv',index=False)
    weights.to_csv(SUMMARY/'selection-weights.csv',index=False)
    summary=selected.groupby(GROUP+['rule']).agg(n=('absolute_error','size'),
        mae=('absolute_error','mean'),median_absolute_error=('absolute_error','median'),
        q25_absolute_error=('absolute_error',lambda x:x.quantile(.25)),
        q75_absolute_error=('absolute_error',lambda x:x.quantile(.75)),
        mean_cv_mse=('cv_mse','mean'),mean_target_mse=('target_mse','mean'),
        mean_gap=('gap','mean'),mean_relative_error=('relative_error','mean')).reset_index()
    summary.to_csv(SUMMARY/'selection-summary.csv',index=False)
    weights.groupby(GROUP+['rule','scheme']).weight.mean().reset_index().to_csv(SUMMARY/'selection-frequencies.csv',index=False)
    candidates=raw.groupby(GROUP+['scheme']).agg(n=('cv_mse','size'),
        mean_cv_mse=('cv_mse','mean'),median_cv_mse=('cv_mse','median'),
        q25_cv_mse=('cv_mse',lambda x:x.quantile(.25)),q75_cv_mse=('cv_mse',lambda x:x.quantile(.75)),
        mean_target_mse=('target_mse','mean'),mean_gap=('gap','mean'),
        mae=('absolute_error','mean'),mean_W_km=('W','mean'),mean_J2=('J2','mean'),
        mean_abs_J2=('J2',lambda x:x.abs().mean()),mean_train=('mean_train','mean'),
        min_train=('min_train','min')).reset_index()
    candidates.to_csv(SUMMARY/'candidate-summary.csv',index=False)
    paired=[]
    for key,frame in selected.groupby(GROUP):
        wide=frame.pivot(index='repeat_id',columns='rule',values='absolute_error')
        delta=wide.W-wide.J2
        paired.append(dict(zip(GROUP,key))|{'W_mae':wide.W.mean(),'J2_mae':wide.J2.mean(),
            'W_minus_J2':delta.mean(),'reduction_percent':100*delta.mean()/wide.W.mean(),
            'delta_q025':delta.quantile(.025),'delta_q975':delta.quantile(.975),
            'J2_better_repetitions':int((delta>1e-12).sum()),
            'W_better_repetitions':int((delta< -1e-12).sum()),
            'equal_repetitions':int((delta.abs()<=1e-12).sum())})
    pd.DataFrame(paired).to_csv(SUMMARY/'paired-comparison.csv',index=False)
    fits['is_full']=fits['fit']=='full'
    fit_summary=fits.groupby(GROUP+['is_full']).agg(n_fits=('n','size'),
        min_signal=('signal','min'),median_signal=('signal','median'),max_signal=('signal','max'),
        min_range_km=('range_km','min'),median_range_km=('range_km','median'),max_range_km=('range_km','max'),
        zero_nugget_count=('zero_nugget','sum'),negligible_signal_count=('negligible_signal','sum'),
        long_range_count=('long_range','sum'),jitter_count=('jitter',lambda x:(x>0).sum()),
        warning_count=('warnings',lambda x:x.notna().sum())).reset_index()
    fit_summary.to_csv(SUMMARY/'fit-summary.csv',index=False)
    # Diagnostic directions are descriptive, not tests of a repeated-field bias.
    signals=[]
    for key,frame in raw.groupby(GROUP):
        for criterion in ['W','J2']:
            x=frame[criterion] if criterion=='W' else frame[criterion].abs()
            signals.append(dict(zip(GROUP,key))|{'diagnostic':criterion,
                'spearman_with_absolute_gap':x.corr(frame.absolute_error,method='spearman')})
    pd.DataFrame(signals).to_csv(SUMMARY/'descriptive-diagnostic-association.csv',index=False)
    # Resample spatial groups of held-out errors, keeping all fitted predictions
    # and CV scores fixed. This is a sensitivity range, not an inferential CI.
    resampled=[]; blocks=[]
    for gi,g in enumerate(x for x in details if x['repeat_id']==1 and x['target']!='uniform'):
        dataset,target_type=g['dataset'],g['target']
        d=pd.read_csv(ROOT/'revision/data/jura-sites.csv')
        site_ids=d.iloc[np.asarray(g['target_rows'])-1].site_id.to_numpy()
        coords=np.asarray(g['target_xy'])
        k=max(4,min(10,len(site_ids)//10))
        labels=KMeans(n_clusters=k,random_state=1400000+gi,n_init=20).fit_predict(coords)
        blocks.extend({'dataset':dataset,'target':target_type,'site_id':site,'block':int(label)}
                      for site,label in zip(site_ids,labels))
        rng=np.random.default_rng(1500000+gi)
        sampled=rng.integers(0,k,size=(2000,k))
        multiplicities=np.asarray([(row[:,None]==np.arange(k)).sum(axis=0)[labels] for row in sampled])
        for model in ['ordinary','trend']:
            p=pred[(pred.dataset==dataset)&(pred.target==target_type)&(pred.repeat_id==1)&
                   (pred.model==model)&(pred.role=='target')].set_index('site_id').loc[site_ids]
            mse=multiplicities@p.squared_error.to_numpy()/multiplicities.sum(axis=1)
            frames=[frame for _,frame in raw[(raw.dataset==dataset)&(raw.target==target_type)&
                                             (raw.model==model)].groupby('repeat_id')]
            vals={}
            for rule in ['W','J2','fixed_random','fixed_knndm']:
                losses=[]
                for frame in frames:
                    w=selection_weights(frame,rule)
                    losses.append(np.abs(mse[:,None]-frame.cv_mse.to_numpy()[None,:])@w)
                vals[rule]=np.mean(losses,axis=0)
            delta=vals['W']-vals['J2']
            for b in range(len(mse)):
                resampled.append({'dataset':dataset,'target':target_type,'model':model,
                    'draw':b+1,'n_blocks':k,'target_mse':mse[b],
                    'W_mae':vals['W'][b],'J2_mae':vals['J2'][b],'W_minus_J2':delta[b],
                    'random_mae':vals['fixed_random'][b],'knndm_mae':vals['fixed_knndm'][b],
                    'J2_minus_random':vals['J2'][b]-vals['fixed_random'][b]})
    pd.DataFrame(blocks).to_csv(SUMMARY/'target-spatial-blocks.csv',index=False)
    boot=pd.DataFrame(resampled)
    boot.to_csv(SUMMARY/'target-block-sensitivity-draws.csv',index=False)
    boot.groupby(GROUP).agg(n_blocks=('n_blocks','first'),
        target_mse_q025=('target_mse',lambda x:x.quantile(.025)),
        target_mse_q975=('target_mse',lambda x:x.quantile(.975)),
        delta_q025=('W_minus_J2',lambda x:x.quantile(.025)),
        delta_q975=('W_minus_J2',lambda x:x.quantile(.975)),
        fraction_J2_better=('W_minus_J2',lambda x:(x>1e-12).mean()),
        J2_minus_random_q025=('J2_minus_random',lambda x:x.quantile(.025)),
        J2_minus_random_q975=('J2_minus_random',lambda x:x.quantile(.975)),
        fraction_random_better=('J2_minus_random',lambda x:(x>1e-12).mean())).reset_index().to_csv(
            SUMMARY/'target-block-sensitivity.csv',index=False)
    print(json.dumps(checks,indent=2))
    print('Wrote Jura verification and target-cluster sensitivity tables.')

if __name__=='__main__':
    main()
