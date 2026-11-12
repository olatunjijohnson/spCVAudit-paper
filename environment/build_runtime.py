#!/usr/bin/env python3
"""Explicitly download and build a pinned R runtime in a separate library.

Requires Python 3.9+, R, a compiler, and the system libraries required by the
locked R packages. Base packages remain supplied by R itself. This does not
modify a user-wide R library or run scientific studies.
"""
import argparse
import csv
import copy
import datetime as dt
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tarfile
import urllib.error
import urllib.request


def utc():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def dcf(raw):
    fields = {}
    key = None
    for line in raw.decode().splitlines():
        if not line.strip():
            continue
        if line[:1].isspace():
            fields[key] += " " + line.strip()
        else:
            key, value = line.split(":", 1)
            fields[key] = value.strip()
    return fields


def dependencies(row):
    result = set()
    for field in ("Depends", "Imports", "LinkingTo"):
        for item in row[field].split(","):
            match = re.match(r"\s*([A-Za-z][A-Za-z0-9.]*)", item)
            if match and match[1] != "R":
                result.add(match[1])
    return result


def source_matches(path, row):
    try:
        with tarfile.open(path, "r:gz") as archive:
            data = dcf(archive.extractfile(row["Package"] + "/DESCRIPTION").read())
        return data["Package"] == row["Package"] and data["Version"] == row["Version"]
    except (OSError, KeyError, tarfile.TarError, AttributeError, UnicodeError):
        return False


def read_source_hashes(path, records):
    if path is None:
        return {}
    rows = list(csv.DictReader(path.open()))
    hashes = {row["Package"]: row for row in rows}
    required = {p for p, row in records.items() if row["Priority"] != "base"}
    if len(hashes) != len(rows) or set(hashes) != required:
        raise ValueError("Source hashes must contain every non-base package exactly once")
    for package, row in hashes.items():
        if row["Version"] != records[package]["Version"] or not re.fullmatch(r"[a-f0-9]{64}", row["sha256"]):
            raise ValueError("Invalid source hash or version: " + package)
    return hashes


def verify_source_hash(archive, package, expected):
    actual = sha256(archive)
    if package in expected and actual != expected[package]["sha256"]:
        raise ValueError("Source archive differs from recorded SHA-256: " + package)
    return actual


def github_pin(provenance):
    commit = provenance.get("RemoteSha") or provenance.get("GithubSHA1")
    if not commit:
        return None
    repository = provenance.get("RemoteUrl", "").removesuffix(".git")
    if not repository and provenance.get("GithubUsername") and provenance.get("GithubRepo"):
        repository = "https://github.com/" + provenance["GithubUsername"] + "/" + provenance["GithubRepo"]
    if not re.fullmatch(r"[a-fA-F0-9]{40}", commit) or not re.fullmatch(
            r"https://github[.]com/[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+", repository):
        raise ValueError("Unsupported recorded source revision: " + str(provenance))
    return repository, commit


def github_source(archive, row, pin, source_cache, cache_only=False):
    repository, commit = pin
    sidecar = archive.with_name(archive.name + ".source.json")
    candidates = [(archive, sidecar)]
    if source_cache:
        candidates.append((source_cache / archive.name, source_cache / sidecar.name))
    for candidate, metadata in candidates:
        if candidate.is_file() and metadata.is_file():
            record = json.loads(metadata.read_text())
            if (record.get("commit") == commit and record.get("repository") == repository
                    and record.get("source_sha256") == sha256(candidate) and source_matches(candidate, row)):
                if candidate.resolve() != archive.resolve():
                    shutil.copy2(candidate, archive)
                    shutil.copy2(metadata, sidecar)
                upstream_name = row["Package"] + "-" + commit + ".github.tar.gz"
                upstream = candidate.parent / upstream_name
                if upstream.is_file():
                    if sha256(upstream) != record["upstream_archive_sha256"]:
                        raise ValueError("Cached upstream Git archive has changed")
                    destination = archive.parent / upstream_name
                    if upstream.resolve() != destination.resolve():
                        shutil.copy2(upstream, destination)
                return record
    if cache_only:
        raise RuntimeError("No verified cached Git source; downloads disabled: " + row["Package"])
    url = "https://codeload.github.com/" + repository.removeprefix("https://github.com/") + "/tar.gz/" + commit
    original = archive.parent / (row["Package"] + "-" + commit + ".github.tar.gz")
    with urllib.request.urlopen(url, timeout=60) as response, original.open("wb") as output:
        shutil.copyfileobj(response, output)
    # Normalize only the archive's outer directory. Keep every source file's
    # bytes unchanged; no package version or implementation is patched.
    with tarfile.open(original, "r:gz") as source, archive.open("wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", filename="", mtime=0) as compressed:
            with tarfile.open(fileobj=compressed, mode="w") as target:
                for member in source.getmembers():
                    parts = Path(member.name).parts
                    if not parts or ".." in parts or Path(member.name).is_absolute() or not (member.isfile() or member.isdir()):
                        raise ValueError("Unexpected upstream archive member: " + member.name)
                    updated = copy.copy(member)
                    updated.name = str(Path(row["Package"], *parts[1:]))
                    updated.pax_headers = {k:v for k,v in member.pax_headers.items() if k != "path"}
                    target.addfile(updated, source.extractfile(member) if member.isfile() else None)
    if not source_matches(archive, row):
        raise ValueError("Recorded Git revision does not contain the locked package version")
    record = {"repository": repository, "commit": commit, "url": url,
              "upstream_archive_sha256": sha256(original), "source_sha256": sha256(archive),
              "preparation": "Outer directory renamed; source-file bytes unchanged"}
    sidecar.write_text(json.dumps(record, indent=2) + "\n")
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("lock", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--source-cache", type=Path)
    parser.add_argument("--provenance", type=Path, help="Recorded repositories/commits for non-CRAN builds")
    parser.add_argument("--source-hashes", type=Path, help="Require the recorded SHA-256 of every non-base source")
    parser.add_argument("--cache-only", action="store_true", help="Disable the builder's downloads; require retained archives and source hashes")
    parser.add_argument("--resume", action="store_true", help="Resume an incomplete build with the identical lock")
    args = parser.parse_args()
    if args.cache_only and args.source_hashes is None:
        parser.error("--cache-only requires --source-hashes")
    lock = args.lock.resolve()
    root = args.destination.resolve()
    rows = list(csv.DictReader(lock.open()))
    records = {row["Package"]: row for row in rows}
    source_hashes = read_source_hashes(args.source_hashes, records)
    provenance = {} if args.provenance is None else {
        row["Package"]: row for row in csv.DictReader(args.provenance.open())}
    for package, row in provenance.items():
        if package not in records or row["Version"] != records[package]["Version"]:
            raise ValueError("Source provenance differs from version lock: " + package)
    if len(records) != len(rows):
        raise ValueError("The lock contains duplicate package names")
    base_versions = {row["Version"] for row in rows if row["Priority"] == "base"}
    actual_R = subprocess.check_output(
        ["Rscript", "--vanilla", "-e", "cat(as.character(getRversion()))"], text=True).strip()
    if base_versions != {actual_R}:
        raise RuntimeError("This lock requires R " + ", ".join(sorted(base_versions)) +
                           "; available R is " + actual_R + ". No packages were built.")
    missing = set().union(*(dependencies(row) for row in rows)) - set(records)
    if missing:
        raise ValueError("Incomplete dependency closure: " + ", ".join(sorted(missing)))
    pending = {p for p, row in records.items() if row["Priority"] != "base"}
    order = []
    while pending:
        ready = sorted(p for p in pending if not (dependencies(records[p]) & pending))
        if not ready:
            raise ValueError("Dependency graph has a cycle")
        order.extend(ready)
        pending.difference_update(ready)

    receipt = root / "BUILD.json"
    if root.exists():
        if not args.resume or not receipt.is_file():
            raise FileExistsError("Choose a new destination or explicitly resume an incomplete build")
        build = json.loads(receipt.read_text())
        if build["lock_sha256"] != sha256(lock) or build["status"] == "complete":
            raise RuntimeError("Cannot resume a different lock or a completed runtime")
    else:
        root.mkdir(parents=True)
        build = {"started_utc": utc(), "lock_sha256": sha256(lock), "packages": {}, "attempts": []}
        shutil.copy2(lock, root / "runtime-lock.csv")
    for directory in ("library", "site-library", "sources", "logs"):
        (root / directory).mkdir(exist_ok=True)
    library = root / "library"
    env = os.environ.copy()
    env.update(R_LIBS_USER=str(library), R_LIBS_SITE=str(root / "site-library"),
               R_ENVIRON_USER="/dev/null", R_PROFILE_USER="/dev/null",
               OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", MAKEFLAGS="-j1")
    env.pop("R_LIBS", None)
    build["attempts"].append({"started_utc": utc(), "script_sha256": sha256(Path(__file__)),
                              "slurm_job": env.get("SLURM_JOB_ID"),
                              "builder_downloads_disabled": args.cache_only,
                              "R_version": subprocess.check_output(["R", "--version"], env=env, text=True).splitlines()[0]})
    snapshot = root / ("builder-attempt" + str(len(build["attempts"])) + ".py")
    if snapshot.exists() and sha256(snapshot) != sha256(Path(__file__)):
        raise RuntimeError("Builder source snapshot already exists with different contents")
    shutil.copy2(Path(__file__), snapshot)
    build["attempts"][-1]["script_snapshot"] = str(snapshot.relative_to(root))
    build["status"] = "running"
    build.pop("error", None)
    build.pop("finished_utc", None)
    if args.provenance:
        provenance_hash = sha256(args.provenance)
        if build.get("source_provenance_sha256", provenance_hash) != provenance_hash:
            raise RuntimeError("Source provenance changed during the isolated build")
        build["source_provenance_sha256"] = provenance_hash
        shutil.copy2(args.provenance, root / "runtime-source-provenance.csv")
    elif build.get("source_provenance_sha256"):
        raise RuntimeError("Resume requires the recorded source-provenance file")
    if args.source_hashes:
        hashes_hash = sha256(args.source_hashes)
        if build.get("source_hashes_sha256", hashes_hash) != hashes_hash:
            raise RuntimeError("Expected source hashes changed during the isolated build")
        build["source_hashes_sha256"] = hashes_hash
        shutil.copy2(args.source_hashes, root / "runtime-source-hashes.csv")
    elif build.get("source_hashes_sha256"):
        raise RuntimeError("Resume requires the recorded source hashes")

    def save():
        receipt.write_text(json.dumps(build, indent=2) + "\n")

    save()
    try:
        for i, package in enumerate(order, 1):
            row = records[package]
            pin = github_pin(provenance.get(package, {}))
            prior = build["packages"].get(package, {})
            description = library / package / "DESCRIPTION"
            if prior.get("exit_code") == 0 and description.is_file():
                installed = dcf(description.read_bytes())
                if installed.get("Version") != row["Version"]:
                    raise RuntimeError("A completed installation has changed: " + package)
                if pin and prior.get("git_commit") != pin[1]:
                    raise RuntimeError("Already-built package needs a different recorded source; use a fresh destination: " + package)
                if package in source_hashes:
                    retained = root / "sources" / (package + "_" + row["Version"] + ".tar.gz")
                    if verify_source_hash(retained, package, source_hashes) != prior["source_sha256"]:
                        raise RuntimeError("Previously built source has changed: " + package)
                print(utc(), "Already built", package, row["Version"], flush=True)
                continue
            filename = package + "_" + row["Version"] + ".tar.gz"
            archive = root / "sources" / filename
            source_url = None
            git_source = github_source(archive, row, pin, args.source_cache, args.cache_only) if pin else None
            if git_source:
                source_url = git_source["url"]
            if not source_matches(archive, row) and args.source_cache:
                cached = args.source_cache / filename
                if source_matches(cached, row):
                    shutil.copy2(cached, archive)
                    source_url = "existing source archive: " + str(cached)
            if not source_matches(archive, row):
                if args.cache_only:
                    raise RuntimeError("No valid cached source; downloads disabled: " + package)
                errors = []
                for url in ("https://cloud.r-project.org/src/contrib/Archive/" + package + "/" + filename,
                            "https://cloud.r-project.org/src/contrib/" + filename):
                    try:
                        with urllib.request.urlopen(url, timeout=60) as response, archive.open("wb") as output:
                            shutil.copyfileobj(response, output)
                        if not source_matches(archive, row):
                            raise ValueError("Downloaded source does not match locked package/version")
                        source_url = url
                        break
                    except (OSError, urllib.error.URLError, ValueError) as error:
                        errors.append(str(error))
                else:
                    raise RuntimeError("Cannot retrieve " + filename + ": " + "; ".join(errors))
            log = root / "logs" / (package + "-attempt" + str(len(build["attempts"])) + ".log")
            command = ["R", "CMD", "INSTALL", "--library=" + str(library), str(archive)]
            item = {"version": row["Version"], "source_sha256": verify_source_hash(archive, package, source_hashes),
                    "source": source_url or "retained in this isolated build", "started_utc": utc(),
                    "command": command, "log": str(log.relative_to(root))}
            if git_source:
                item.update(git_commit=pin[1], source_preparation=git_source)
            build["packages"][package] = item
            save()
            print(utc(), "Building", i, "of", len(order), package, row["Version"], flush=True)
            with log.open("w") as output:
                result = subprocess.run(command, env=env, stdout=output, stderr=subprocess.STDOUT)
            item.update(finished_utc=utc(), exit_code=result.returncode, log_sha256=sha256(log))
            save()
            if result.returncode:
                raise RuntimeError("Installation failed: " + package + "; see " + str(log))

        validation = r'''
args <- commandArgs(trailingOnly = TRUE)
lock <- read.csv(args[1], stringsAsFactors = FALSE)
library <- normalizePath(args[2])
.libPaths(c(library, .Library))
stopifnot(all(vapply(lock$Package, requireNamespace, logical(1), quietly = TRUE)))
paths <- vapply(lock$Package, find.package, "")
actual <- vapply(lock$Package, function(p) as.character(packageVersion(p)), "")
stopifnot(all(package_version(actual) == package_version(lock$Version)))
nonbase <- lock$Priority != "base"
stopifnot(all(dirname(normalizePath(paths[nonbase])) == library))
write.csv(data.frame(package = lock$Package, version = actual, path = paths),
          file.path(args[3], "resolved-packages.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(args[3], "runtime-session.txt"))
cat("Every non-base dependency resolves from the isolated library.\n")
'''
        validation_script = root / "validate.R"
        validation_script.write_text(validation)
        with (root / "logs/validation.log").open("w") as output:
            subprocess.run(["Rscript", "--vanilla", str(validation_script), str(lock), str(library), str(root)],
                           env=env, stdout=output, stderr=subprocess.STDOUT, check=True)
        build.update(status="complete", finished_utc=utc(), nonbase_packages_built=len(order),
                     base_packages_supplied_by_R=len(records)-len(order))
        build["attempts"][-1].update(status="complete", finished_utc=build["finished_utc"])
    except Exception as error:
        build.update(status="failed", finished_utc=utc(), error=str(error))
        build["attempts"][-1].update(status="failed", finished_utc=build["finished_utc"], error=str(error))
        save()
        raise
    save()
    print("Completed isolated dependency build:", receipt, flush=True)


if __name__ == "__main__":
    main()
