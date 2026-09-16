#!/usr/bin/env python3
"""
Harvest spatial transcriptomics and single-cell RNA-seq dataset records from
NCBI GEO (E-utilities) and EBI ArrayExpress (BioStudies API), classify them,
and write data/datasets.json + data/meta.json.

Only stdlib is used so the GitHub Action needs no dependency install.

Provenance note
---------------
Accession, title, description, organism, sample count, release date and PubMed
ID come verbatim from the source databases. Everything else (modality,
platform, platform class, disease category, tissue, cell count) is DERIVED by
the rules in this file from the title + description text, and is therefore a
best-effort label, not a curator's judgement. The site labels these fields as
auto-classified.
"""

import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
BIOSTUDIES = "https://www.ebi.ac.uk/biostudies/api/v1/arrayexpress/search"
NCBI_KEY = os.environ.get("NCBI_API_KEY", "").strip()
UA = "spatial-sc-atlas/1.0 (https://github.com/cemalyildiz/spatial-sc-atlas)"

# Polite pacing: 10 req/s with an API key, 3 req/s without.
PAUSE = 0.12 if NCBI_KEY else 0.40


# --------------------------------------------------------------------------
# HTTP helpers
# --------------------------------------------------------------------------
def _request(url, data=None, tries=4):
    last = None
    for attempt in range(tries):
        try:
            req = urllib.request.Request(
                url,
                data=data.encode() if isinstance(data, str) else data,
                headers={
                    "User-Agent": UA,
                    "Accept": "application/json",
                    **(
                        {"Content-Type": "application/x-www-form-urlencoded"}
                        if data
                        else {}
                    ),
                },
            )
            with urllib.request.urlopen(req, timeout=90) as resp:
                return json.loads(resp.read().decode("utf-8", "replace"))
        except Exception as exc:  # noqa: BLE001 - network flakiness is expected
            last = exc
            time.sleep(2 * (attempt + 1))
    print(f"  ! giving up on {url[:90]}: {last}", file=sys.stderr)
    return None


def esearch(term, retmax=6000):
    params = {
        "db": "gds",
        "retmode": "json",
        "retmax": str(retmax),
        "term": term,
    }
    if NCBI_KEY:
        params["api_key"] = NCBI_KEY
    j = _request(EUTILS + "esearch.fcgi?" + urllib.parse.urlencode(params))
    time.sleep(PAUSE)
    if not j:
        return []
    return j.get("esearchresult", {}).get("idlist", []) or []


def esummary(uids, batch=250):
    out = []
    for i in range(0, len(uids), batch):
        chunk = uids[i : i + batch]
        body = {"db": "gds", "retmode": "json", "id": ",".join(chunk)}
        if NCBI_KEY:
            body["api_key"] = NCBI_KEY
        j = _request(EUTILS + "esummary.fcgi", data=urllib.parse.urlencode(body))
        time.sleep(PAUSE)
        if not j or "result" not in j:
            continue
        res = j["result"]
        for uid in res.get("uids", []):
            d = res.get(uid)
            if not d:
                continue
            out.append(
                {
                    "acc": d.get("accession", ""),
                    "title": d.get("title", "") or "",
                    "summary": d.get("summary", "") or "",
                    "taxon": d.get("taxon", "") or "",
                    "gdstype": d.get("gdstype", "") or "",
                    "n": d.get("n_samples", 0) or 0,
                    "pdat": d.get("pdat", "") or "",
                    "pmid": (d.get("pubmedids") or [""])[0],
                    "ftp": d.get("ftplink", "") or "",
                }
            )
        print(f"  esummary {min(i + batch, len(uids))}/{len(uids)}", flush=True)
    return out


def ae_search(query, max_pages=30):
    out = {}
    for page in range(1, max_pages + 1):
        url = (
            BIOSTUDIES
            + "?"
            + urllib.parse.urlencode(
                {"query": query, "pageSize": "100", "page": str(page)}
            )
        )
        j = _request(url)
        time.sleep(0.25)
        if not j:
            break
        hits = j.get("hits") or []
        if not hits:
            break
        for h in hits:
            acc = h.get("accession", "")
            if acc:
                out[acc] = {
                    "acc": acc,
                    "title": h.get("title", "") or "",
                    "content": (h.get("content", "") or "")[:900],
                    "rdate": h.get("release_date", "") or "",
                    "files": h.get("files", 0) or 0,
                }
        print(f"  arrayexpress page {page}: {len(out)} total", flush=True)
        if len(hits) < 100:
            break
    return list(out.values())


# --------------------------------------------------------------------------
# Classification vocabularies
# --------------------------------------------------------------------------
def _c(patterns):
    return [(name, re.compile(rx, re.I)) for name, rx in patterns]


IMAGING = _c(
    [
        ("Xenium", r"\bxenium\b"),
        ("MERFISH / MERSCOPE", r"\bmerfish\b|\bmerscope\b|\bvizgen\b"),
        ("CosMx", r"\bcosmx\b|spatial molecular imager"),
        ("seqFISH", r"\bseqfish\b"),
        ("STARmap", r"\bstarmap\b"),
        ("In situ sequencing", r"\bin situ sequencing\b|\bhybiss\b|\bcartana\b"),
        ("IMC / MIBI", r"\bmibi-?tof\b|imaging mass cytometry"),
        ("Molecular Cartography", r"molecular cartography|resolve bioscience"),
    ]
)

SEQUENCING = _c(
    [
        ("Visium HD", r"visium\s*hd"),
        ("Visium", r"\bvisium\b"),
        ("GeoMx DSP", r"\bgeomx\b|digital spatial profil"),
        ("Slide-seq", r"\bslide-?seq\b"),
        ("Stereo-seq", r"\bstereo-?seq\b|\bstereoseq\b"),
        ("Slide-tags", r"\bslide-?tags\b"),
        ("DBiT-seq", r"\bdbit-?seq\b"),
        ("Curio Seeker", r"\bcurio seeker\b"),
        ("Open-ST", r"\bopen-?st\b"),
        (
            "Spatial (unspecified)",
            r"spatial transcriptom|spatially resolved|spatial gene expression|spatial rna",
        ),
    ]
)

SC_RE = re.compile(
    r"single[- ]cell rna|scrna-?seq|snrna-?seq|single[- ]nucleus rna|"
    r"single[- ]cell transcriptom|10x genomics chromium|\bchromium\b|cite-?seq|"
    r"single[- ]cell sequencing|smart-?seq2?|single[- ]cell atlas",
    re.I,
)

SC_PLATFORM = _c(
    [
        ("10x Chromium", r"\bchromium\b|10x genomics"),
        ("Smart-seq2", r"smart-?seq2?"),
        ("snRNA-seq", r"snrna-?seq|single[- ]nucleus"),
        ("BD Rhapsody", r"bd rhapsody"),
        ("Drop-seq", r"drop-?seq"),
        ("CITE-seq", r"cite-?seq"),
    ]
)

# Cancer subtypes. A record can carry several (LUAD *and* EGFR-mutant), so this
# is matched as a multi-value field rather than a single winner.
#
# Every pattern here must be unambiguous IN THIS CORPUS. Short abbreviations are
# the trap: ILC is "innate lymphoid cell" far more often than "invasive lobular
# carcinoma", FL and MCL collide with gene and cell-type names, IDC is also
# "interdigitating dendritic cell", BCC is also "B-cell ...". Those are matched
# by their full phrase only. Add a bare abbreviation here only after checking it
# against the catalog for collisions.
SUBTYPES = _c(
    [
        # lung
        ("LUAD", r"\bluad\b|lung adenocarcinoma|adenocarcinoma of (the )?lungs?|pulmonary adenocarcinoma"),
        ("LUSC", r"\blusc\b|lung squamous|squamous cell (carcinoma )?of (the )?lungs?|pulmonary squamous"),
        ("LCNEC / large cell", r"\blcnec\b|large[- ]cell (neuroendocrine )?carcinoma"),
        ("Adenosquamous", r"adenosquamous"),
        ("EGFR-mutant", r"\begfr[- ](mutant|mutated|driven|altered)|egfr mutation"),
        ("ALK-rearranged", r"\balk[- ](rearrang|fusion|positive|translocat)"),
        ("KRAS-mutant", r"\bkras[- ]?(g12|g13|mutant|mutated|mutation)"),
        # breast
        ("TNBC", r"\btnbc\b|triple[- ]negative breast"),
        ("HER2-positive", r"her2[- ]?(positive|enriched|amplified)|erbb2[- ]amplif"),
        ("Luminal", r"luminal [ab]\b|luminal subtype|luminal breast"),
        ("Invasive ductal (IDC)", r"invasive ductal carcinoma"),
        ("Invasive lobular (ILC)", r"invasive lobular carcinoma"),
        ("DCIS", r"\bdcis\b|ductal carcinoma in situ"),
        # kidney
        ("ccRCC", r"\bccrcc\b|clear[- ]cell renal"),
        ("pRCC", r"\bprcc\b|papillary renal"),
        ("chRCC", r"\bchrcc\b|chromophobe renal"),
        # prostate
        ("CRPC", r"\bm?crpc\b|castration[- ]resistant"),
        ("NEPC", r"\bnepc\b|neuroendocrine prostate"),
        # brain
        ("IDH-mutant", r"\bidh[12]?[- ]?(mutant|mutated)"),
        ("IDH-wildtype", r"\bidh[12]?[- ]?(wild[- ]?type|wt)\b"),
        ("DIPG / diffuse midline", r"\bdipg\b|diffuse (intrinsic pontine|midline) glioma"),
        # blood
        ("AML", r"\baml\b|acute myeloid leuk"),
        ("B-ALL", r"\bb-all\b|b[- ]cell acute lymphoblastic"),
        ("T-ALL", r"\bt-all\b|t[- ]cell acute lymphoblastic"),
        ("CLL", r"\bcll\b|chronic lymphocytic leuk"),
        ("CML", r"\bcml\b|chronic myeloid leuk"),
        ("MDS", r"\bmds\b|myelodysplastic"),
        ("DLBCL", r"\bdlbcl\b|diffuse large b[- ]cell lymphoma"),
        ("Follicular lymphoma", r"follicular lymphoma"),
        ("Mantle cell lymphoma", r"mantle cell lymphoma"),
        ("Hodgkin", r"hodgkin"),
        # GI
        ("ESCC", r"\bescc\b|(o)?esophageal squamous"),
        ("Esophageal adenocarcinoma", r"(o)?esophageal adenocarcinoma"),
        ("iCCA", r"\bicca\b|intrahepatic cholangiocarcinoma"),
        ("eCCA", r"\becca\b|extrahepatic cholangiocarcinoma|perihilar cholangio"),
        ("MSI-high", r"\bmsi-?h\b|microsatellite instab"),
        ("MSS", r"microsatellite stable"),
        ("Signet ring", r"signet[- ]ring"),
        # pancreas / liver
        ("PDAC", r"\bpdac\b|pancreatic ductal adenocarcinoma"),
        ("HCC", r"\bhcc\b|hepatocellular carcinoma"),
        # bladder
        ("MIBC", r"\bmibc\b|muscle[- ]invasive bladder"),
        ("NMIBC", r"\bnmibc\b|non[- ]muscle[- ]invasive bladder"),
        # gynaecological
        ("HGSOC", r"\bhgsoc\b|high[- ]grade serous"),
        ("Ovarian clear cell", r"ovarian clear cell"),
        # head & neck / skin
        ("HNSCC", r"\bhnscc\b|head and neck squamous"),
        ("OSCC", r"\boscc\b|oral squamous cell"),
        ("HPV-positive", r"hpv[- ]?(positive|associated|driven)"),
        ("Basal cell carcinoma", r"basal cell carcinoma"),
        ("cSCC", r"\bcscc\b|cutaneous squamous cell"),
    ]
)

CANCER = _c(
    [
        ("NSCLC", r"non-?small[- ]cell lung|\bnsclc\b|lung adenocarcinoma|\bluad\b|lung squamous|\blusc\b"
                  r"|adenocarcinoma of (the )?lungs?|pulmonary adenocarcinoma|\blcnec\b|bronchioloalveolar"
                  r"|lepidic|adenosquamous carcinoma of (the )?lung|large[- ]cell carcinoma of (the )?lung"),
        ("SCLC", r"small[- ]cell lung cancer|\bsclc\b"),
        ("Lung cancer (other)", r"lung (cancer|tumou?r|carcinoma|metasta|neoplas)|pulmonary (cancer|carcinoma)|mesothelioma"),
        ("Breast cancer", r"breast (cancer|tumou?r|carcinoma)|triple[- ]negative breast|\btnbc\b"
                          r"|ductal carcinoma in situ|\bdcis\b|invasive (ductal|lobular) carcinoma"
                          r"|luminal [ab] (breast|subtype)"),
        ("Colorectal cancer", r"colorectal|\bcrc\b|colon (cancer|adenocarcinoma|tumou?r)|rectal cancer"),
        ("Glioma / GBM", r"glioblastoma|\bgbm\b|\bglioma\b|astrocytoma|oligodendroglioma|medulloblastoma"
                         r"|ependymoma|\bdipg\b|diffuse (intrinsic pontine|midline) glioma"),
        ("Meningioma", r"meningioma"),
        ("Pancreatic cancer", r"pancreatic (cancer|ductal|adenocarcinoma|tumou?r)|\bpdac\b"),
        ("Prostate cancer", r"prostate (cancer|tumou?r|adenocarcinoma)|\bm?crpc\b|\bnepc\b|castration[- ]resistant"),
        ("Melanoma", r"melanoma"),
        ("Liver / HCC", r"hepatocellular|\bhcc\b|liver (cancer|tumou?r)|cholangiocarcinoma|\bicca\b|\becca\b"),
        ("Gastric cancer", r"gastric (cancer|adenocarcinoma|tumou?r)|stomach cancer"),
        ("Ovarian cancer", r"ovarian (cancer|carcinoma|tumou?r|clear cell)|\bhgsoc\b"),
        ("Kidney / RCC", r"renal cell carcinoma|\brcc\b|\bccrcc\b|\bprcc\b|\bchrcc\b|kidney (cancer|tumou?r)|wilms"),
        ("Head & neck cancer", r"head and neck|\bhnscc\b|oral squamous|\boscc\b|nasopharyngeal|laryngeal cancer"),
        ("Lymphoma", r"lymphoma|\bdlbcl\b"),
        ("Leukemia", r"leukemia|leukaemia|\baml\b|\bcll\b|\bcml\b|\bapl\b|\bb-all\b|\bt-all\b|myelodysplas"),
        ("Multiple myeloma", r"multiple myeloma|plasma cell myeloma"),
        ("Esophageal cancer", r"(o)?esophageal (cancer|carcinoma|adenocarcinoma|squamous)|\bescc\b"),
        ("Bladder cancer", r"bladder cancer|urothelial (carcinoma|cancer)|\bmibc\b|\bnmibc\b"),
        ("Cervical cancer", r"cervical (cancer|carcinoma)"),
        ("Endometrial cancer", r"endometrial (cancer|carcinoma)|uterine (cancer|carcinoma)"),
        ("Sarcoma", r"sarcoma|\bgist\b"),
        ("Neuroblastoma", r"neuroblastoma"),
        ("Thyroid cancer", r"thyroid (cancer|carcinoma)"),
        ("Skin cancer (non-melanoma)", r"basal cell carcinoma|cutaneous squamous|\bcscc\b"),
        ("Cancer (other / pan-cancer)", r"\bcancers?\b|\btumou?rs?\b|carcinoma|malignan|metasta|neoplas|oncogen|oncolog"),
    ]
)

OTHER_DISEASE = _c(
    [
        ("Neurodegenerative", r"alzheimer|parkinson|\bals\b|amyotrophic lateral|huntington|multiple sclerosis|dementia|neurodegener|tauopath"),
        ("Fibrosis", r"\bfibrosis\b|fibrotic|\bipf\b|idiopathic pulmonary fibrosis|cirrhosis"),
        ("Autoimmune / inflammatory", r"rheumatoid|\blupus\b|psoriasis|crohn|ulcerative colitis|inflammatory bowel|\bibd\b|autoimmun|atopic dermatitis|\bsle\b"),
        ("Cardiovascular", r"myocardial infarction|heart failure|cardiomyopath|atheroscleros|cardiac (injury|disease|remodel)|aortic aneurysm"),
        ("Infectious disease", r"covid|sars-cov-2|tuberculosis|influenza|\bhiv\b|\bsepsis\b|malaria|helicobacter|viral infection"),
        ("Metabolic", r"diabet|obesity|\bnafld\b|\bnash\b|steatohepatitis|metabolic syndrome|insulin resistance"),
        ("Kidney disease", r"chronic kidney disease|renal fibrosis|nephropath|glomerulo"),
        ("Respiratory (non-cancer)", r"\bcopd\b|\basthma\b|cystic fibrosis|\bards\b|emphysema|bronchopulmonary"),
        ("Injury / regeneration", r"\bwound\b|regenerat|\binjury\b|\bischemi"),
        ("Development / embryo", r"\bembryo|organogenes|\bfetal\b|\bfoetal\b|gestation|developmental atlas"),
        ("Ageing", r"\bageing\b|\baging\b|senescen"),
    ]
)

HEALTHY = re.compile(
    r"healthy (donor|tissue|control|human|adult)|normal (tissue|human|donor)|"
    r"reference (atlas|map)|cell atlas",
    re.I,
)

TISSUE = _c(
    [
        ("Lung", r"\blungs?\b|pulmonary|bronch|alveol|\bairway"),
        ("Brain / CNS", r"\bbrain\b|cortex|cortical|hippocamp|cerebell|\bneuron|\bcns\b|spinal cord|striatum|\bglia"),
        ("Breast", r"\bbreast\b|mammary"),
        ("Colon / intestine", r"\bcolon\b|colorect|intestin|\bileum\b|\brectal\b|duoden|\bcolitis\b"),
        ("Pancreas", r"pancrea|\bislet"),
        ("Liver", r"\bliver\b|hepatic|hepatocyt|\bhcc\b"),
        ("Kidney", r"\bkidney\b|\brenal\b|nephron|glomerul"),
        ("Skin", r"\bskin\b|epiderm|dermal|cutaneous|keratinocyt"),
        ("Prostate", r"prostat"),
        ("Stomach", r"\bstomach\b|\bgastric\b"),
        ("Ovary / uterus / placenta", r"\bovar|uter|endometri|\bcervix\b|\bcervical\b|placenta|\bdecidua"),
        ("Lymph node / spleen", r"lymph node|\bspleen\b|\btonsil|lymphoid"),
        ("Bone marrow / blood", r"bone marrow|\bblood\b|\bpbmc\b|h(a)?ematopoie|\bmyeloid\b"),
        ("Heart", r"\bheart\b|cardiac|myocard|ventric|cardiomyocyt"),
        ("Bladder / urinary", r"\bbladder\b|urothel|\burinary\b"),
        ("Muscle", r"skeletal muscle|\bmyofib|\bmuscles?\b"),
        ("Eye / retina", r"\bretina|\beyes?\b|\bcornea"),
        ("Bone / cartilage", r"\bbones?\b|cartilage|osteo|chondrocyt"),
        ("Testis", r"\btest(is|es|icular)\b|\bsperm"),
        ("Thymus", r"thymus|thymic"),
        ("Adrenal", r"adrenal"),
        ("Thyroid", r"thyroid"),
        ("Esophagus", r"esophag|oesophag"),
        ("Head & neck", r"\boral\b|\btongue\b|salivary|nasopharyn|laryn|\bhnscc\b"),
        ("Adipose", r"adipose"),
        ("Embryo / fetal", r"\bembryo|\bfetal\b|\bfoetal\b"),
    ]
)

# Phrases that contain a cancer keyword but say nothing about the study being
# about cancer. They are blanked out before disease matching, otherwise every
# TNF / TNFR immunology study lands in the pan-cancer bucket.
DISEASE_NOISE = re.compile(
    r"tumou?r necrosis factor|\btnf-?(alpha|a|r\d?)\b|\btnfrsf\d*\b|"
    r"tumou?r necrosis|anti-?tnf",
    re.I,
)

CELLS_RE = re.compile(
    r"([\d]{1,3}(?:[,\.]\d{3})+|\d{4,9})\s*(?:single[- ])?(?:cells|nuclei)\b", re.I
)
CONTROLLED_RE = re.compile(r"\bdbgap\b|\bega\b|controlled access|restricted access", re.I)


def first_match(vocab, text):
    for name, rx in vocab:
        if rx.search(text):
            return name
    return None


def all_matches(vocab, text, limit=6):
    """Every label whose pattern hits — for fields a record can hold several of."""
    out = []
    for name, rx in vocab:
        if rx.search(text):
            out.append(name)
            if len(out) >= limit:
                break
    return out


def best_tissue(text):
    """Pick the tissue with the most keyword hits, not merely the first listed."""
    best, best_score = None, 0
    for name, rx in TISSUE:
        score = len(rx.findall(text))
        if score > best_score:
            best, best_score = name, score
    return best or "Other"


def cell_count(text):
    m = CELLS_RE.search(text)
    if not m:
        return 0
    try:
        v = int(m.group(1).replace(",", "").replace(".", ""))
    except ValueError:
        return 0
    return v if 500 <= v <= 50_000_000 else 0


def organism(taxon_text, fallback_text=""):
    t = taxon_text or fallback_text
    human = bool(re.search(r"Homo sapiens", t, re.I))
    mouse = bool(re.search(r"Mus musculus", t, re.I))
    if human and mouse:
        return "Human & Mouse"
    if human:
        return "Human"
    if mouse:
        return "Mouse"
    return (taxon_text or "Other")[:40] or "Other"


def squash(s, n):
    return re.sub(r"\s+", " ", s or "").strip()[:n]


def classify(title, body, extra=""):
    """Return the derived fields for one record."""
    text = f"{title} || {body} || {extra}"
    img = first_match(IMAGING, text)
    seq = None if img else first_match(SEQUENCING, text)
    platform = img or seq
    has_sc = bool(SC_RE.search(text))
    disease_text = DISEASE_NOISE.sub(" ", text)
    cancer = first_match(CANCER, disease_text)
    other = None if cancer else first_match(OTHER_DISEASE, disease_text)
    if not cancer and not other and HEALTHY.search(text):
        other = "Healthy reference atlas"
    return {
        "platform": platform,
        "platform_class": "imaging" if img else ("sequencing" if seq else None),
        "sc_platform": first_match(SC_PLATFORM, text) if has_sc else None,
        "has_sc": has_sc,
        "disease": cancer or other or "Unspecified",
        "disease_group": "Cancer" if cancer else ("Other disease" if other else "Unspecified"),
        # Subtypes are only meaningful once a record reads as cancer; leaving the
        # gate off would tag immunology studies with blood-cancer abbreviations.
        "subtypes": all_matches(SUBTYPES, disease_text) if cancer else [],
        "tissue": best_tissue(text),
        "cells": cell_count(text),
        "access": "Controlled" if CONTROLLED_RE.search(text) else "Open",
    }


# --------------------------------------------------------------------------
# Queries
# --------------------------------------------------------------------------
SPATIAL_QUERY = (
    '("spatial transcriptomics"[All Fields] OR "spatially resolved transcriptomics"[All Fields] '
    'OR Visium[All Fields] OR Xenium[All Fields] OR MERFISH[All Fields] OR MERSCOPE[All Fields] '
    'OR CosMx[All Fields] OR GeoMx[All Fields] OR "Slide-seq"[All Fields] OR "Stereo-seq"[All Fields] '
    'OR seqFISH[All Fields] OR STARmap[All Fields] OR "DBiT-seq"[All Fields] OR "Slide-tags"[All Fields] '
    'OR "spatial gene expression"[All Fields]) AND gse[Entry Type] '
    'AND ("Homo sapiens"[Organism] OR "Mus musculus"[Organism])'
)

SC_CANCER_QUERY = (
    '("single cell RNA"[All Fields] OR "scRNA-seq"[All Fields] '
    'OR "single-cell RNA sequencing"[All Fields] OR "snRNA-seq"[All Fields]) '
    'AND (cancer[All Fields] OR tumor[All Fields] OR tumour[All Fields] '
    'OR carcinoma[All Fields] OR malignant[All Fields]) AND gse[Entry Type] '
    'AND "Homo sapiens"[Organism]'
)

AE_SPATIAL_QUERY = (
    'spatial transcriptomics OR Visium OR Xenium OR MERFISH OR CosMx OR GeoMx '
    'OR "Slide-seq" OR "Stereo-seq"'
)
AE_SC_QUERY = 'single cell RNA-seq AND (cancer OR tumour OR carcinoma)'


def geo_url(acc):
    return f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}"


def ae_url(acc):
    return f"https://www.ebi.ac.uk/biostudies/arrayexpress/studies/{acc}"


def build():
    rows = []
    seen = set()

    print("GEO: spatial search", flush=True)
    spatial_ids = esearch(SPATIAL_QUERY)
    print(f"  {len(spatial_ids)} ids", flush=True)
    spatial = esummary(spatial_ids)

    for r in spatial:
        c = classify(r["title"], r["summary"], r["gdstype"])
        if not c["platform"]:
            continue
        if r["acc"] in seen:
            continue
        seen.add(r["acc"])
        rows.append(
            {
                "accession": r["acc"],
                "db": "GEO",
                "title": squash(r["title"], 190),
                "description": squash(r["summary"], 320),
                "modality": "paired" if c["has_sc"] else "spatial",
                "platform": c["platform"],
                "platform_class": c["platform_class"],
                "sc_platform": c["sc_platform"],
                "disease": c["disease"],
                "disease_group": c["disease_group"],
                "subtypes": c["subtypes"],
                "tissue": c["tissue"],
                "organism": organism(r["taxon"]),
                "samples": r["n"],
                "cells": c["cells"],
                "year": (r["pdat"] or "")[:4],
                "date": r["pdat"],
                "pmid": str(r["pmid"] or ""),
                "access": c["access"],
                "url": geo_url(r["acc"]),
                "ftp": r.get("ftp", ""),
            }
        )

    print("GEO: single-cell cancer search", flush=True)
    sc_ids = [i for i in esearch(SC_CANCER_QUERY) if i not in set(spatial_ids)]
    print(f"  {len(sc_ids)} new ids", flush=True)
    sc = esummary(sc_ids)

    for r in sc:
        c = classify(r["title"], r["summary"], r["gdstype"])
        if not c["has_sc"] or c["disease_group"] != "Cancer":
            continue
        if r["acc"] in seen:
            continue
        seen.add(r["acc"])
        rows.append(
            {
                "accession": r["acc"],
                "db": "GEO",
                "title": squash(r["title"], 190),
                "description": squash(r["summary"], 320),
                "modality": "singlecell",
                "platform": c["sc_platform"] or "scRNA-seq (unspecified)",
                "platform_class": "sequencing",
                "sc_platform": c["sc_platform"],
                "disease": c["disease"],
                "disease_group": c["disease_group"],
                "subtypes": c["subtypes"],
                "tissue": c["tissue"],
                "organism": organism(r["taxon"]),
                "samples": r["n"],
                "cells": c["cells"],
                "year": (r["pdat"] or "")[:4],
                "date": r["pdat"],
                "pmid": str(r["pmid"] or ""),
                "access": c["access"],
                "url": geo_url(r["acc"]),
                "ftp": r.get("ftp", ""),
            }
        )

    print("ArrayExpress", flush=True)
    ae = {r["acc"]: r for r in ae_search(AE_SPATIAL_QUERY)}
    for r in ae_search(AE_SC_QUERY):
        ae.setdefault(r["acc"], r)

    for r in ae.values():
        # E-GEOD-* records are ArrayExpress mirrors of GEO series; skip duplicates.
        if r["acc"].upper().startswith("E-GEOD"):
            continue
        text = f"{r['title']} {r['content']}"
        if not re.search(r"Homo sapiens|Mus musculus", text, re.I):
            continue
        c = classify(r["title"], r["content"])
        if not c["platform"] and not c["has_sc"]:
            continue
        if not c["platform"] and c["disease_group"] != "Cancer":
            continue
        if r["acc"] in seen:
            continue
        seen.add(r["acc"])
        rows.append(
            {
                "accession": r["acc"],
                "db": "ArrayExpress",
                "title": squash(r["title"], 190),
                "description": squash(r["content"], 320),
                "modality": (
                    ("paired" if c["has_sc"] else "spatial")
                    if c["platform"]
                    else "singlecell"
                ),
                "platform": c["platform"] or c["sc_platform"] or "scRNA-seq (unspecified)",
                "platform_class": c["platform_class"] or "sequencing",
                "sc_platform": c["sc_platform"],
                "disease": c["disease"],
                "disease_group": c["disease_group"],
                "subtypes": c["subtypes"],
                "tissue": c["tissue"],
                "organism": organism("", text),
                "samples": 0,
                "cells": c["cells"],
                "year": (r["rdate"] or "")[:4],
                "date": r["rdate"],
                "pmid": "",
                "access": c["access"],
                "url": ae_url(r["acc"]),
                "ftp": "",
            }
        )

    rows.sort(key=lambda x: (x["date"] or ""), reverse=True)
    return rows


def summarise(rows):
    def tally(key):
        out = {}
        for r in rows:
            out[r[key] or "Unspecified"] = out.get(r[key] or "Unspecified", 0) + 1
        return dict(sorted(out.items(), key=lambda kv: -kv[1]))

    return {
        "generated": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
        "total": len(rows),
        "by_modality": tally("modality"),
        "by_platform_class": tally("platform_class"),
        "by_disease_group": tally("disease_group"),
        "by_database": tally("db"),
        "by_year": tally("year"),
        "sources": [
            "NCBI GEO (E-utilities esearch/esummary, db=gds)",
            "EBI ArrayExpress (BioStudies search API)",
        ],
    }


def main():
    DATA.mkdir(parents=True, exist_ok=True)
    rows = build()

    if len(rows) < 200:
        print(
            f"Harvest returned only {len(rows)} records - refusing to overwrite "
            "the existing catalog with a likely-truncated result.",
            file=sys.stderr,
        )
        return 1

    previous = set()
    existing = DATA / "datasets.json"
    if existing.exists():
        try:
            previous = {
                r["accession"] for r in json.loads(existing.read_text())["datasets"]
            }
        except Exception:  # noqa: BLE001
            previous = set()

    new_accessions = [r["accession"] for r in rows if r["accession"] not in previous]

    meta = summarise(rows)
    meta["new_since_last_run"] = len(new_accessions) if previous else 0

    existing.write_text(
        json.dumps({"meta": meta, "datasets": rows}, ensure_ascii=False, separators=(",", ":"))
    )
    (DATA / "meta.json").write_text(json.dumps(meta, ensure_ascii=False, indent=2))

    print(f"\nWrote {len(rows)} datasets ({meta['new_since_last_run']} new).")
    print(json.dumps(meta["by_modality"], indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
