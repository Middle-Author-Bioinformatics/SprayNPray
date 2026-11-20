#!/usr/bin/env python3
"""Create an offline, portable assembly/ANI/phylogeny dashboard from hybrid_results.

Python 3.9+ standard library only. Generates dashboard.html, downloadable
assemblies and summary tables, and a ZIP bundle. Does not run analyses.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import io
import json
import math
import re
import shutil
import zipfile
from datetime import datetime, timezone
from pathlib import Path

VERSION='1.0.0'
FASTA_SUFFIXES=('.fa','.fna','.fasta','.fa.gz','.fna.gz','.fasta.gz')


def read_table(path,delimiter=','):
    if not path.is_file():
        return []
    csv.field_size_limit(2**31-1)
    with path.open(newline='',encoding='utf-8-sig') as f:
        return list(csv.DictReader(f,delimiter=delimiter))


def numeric(value):
    try:
        x=float(str(value).replace(',',''))
        return x if math.isfinite(x) else None
    except (TypeError,ValueError):
        return None


def norm(label):
    return re.sub('_+','_',label)


def fasta_stats(path):
    lengths=[];gc=acgt=ambiguous=ns=0;header=None;length=0;circular=0
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt') as f:
        for line in f:
            line=line.strip()
            if not line: continue
            if line.startswith('>'):
                if header is not None:
                    if not length: raise ValueError('Empty FASTA record in '+str(path))
                    lengths.append(length)
                header=line[1:];length=0
                circular+=bool(re.search(r'\bcircular[=:](?:true|yes|1)\b',header,re.I))
            else:
                if header is None: raise ValueError('Sequence before FASTA header: '+str(path))
                seq=''.join(line.split()).upper()
                if not re.fullmatch('[ACGTRYSWKMBDHVN]+',seq): raise ValueError('Non-DNA FASTA: '+str(path))
                length+=len(seq);gc+=seq.count('G')+seq.count('C')
                known=sum(seq.count(x) for x in 'ACGT');acgt+=known;ambiguous+=len(seq)-known;ns+=seq.count('N')
        if header is not None:
            if not length: raise ValueError('Empty FASTA record in '+str(path))
            lengths.append(length)
    if not lengths: raise ValueError('No FASTA records: '+str(path))
    total=sum(lengths);acc=0;n50=0
    for length in sorted(lengths,reverse=True):
        acc+=length
        if acc>=total/2: n50=length;break
    return dict(length_bp=total,contigs=len(lengths),n50=n50,largest_contig=max(lengths),
                gc_percent=100*gc/acgt if acgt else None,ambiguous_bp=ambiguous,
                n_per_100kb=100000*ns/total,circular_contigs=circular)


def quast_metrics(path,sample):
    if not path.is_file(): return {}
    with path.open() as f: rows=list(csv.reader(f,delimiter='\t'))
    if not rows or len(rows[0])<2: return {}
    columns=rows[0][1:]
    index=columns.index(sample)+1 if sample in columns else 1 if len(columns)==1 else None
    return {r[0]:r[index] for r in rows[1:] if index is not None and len(r)>index}


def sample_files(root):
    combined={}
    folder=root/'assemblies_combined'
    if folder.is_dir():
        for p in sorted(folder.iterdir()):
            suffix=next((s for s in sorted(FASTA_SUFFIXES,key=len,reverse=True) if p.name.endswith(s)),None)
            if p.is_file() and suffix:
                sample=p.name[:-len(suffix)]
                if sample in combined: raise ValueError('Duplicate combined assembly: '+sample)
                combined[sample]=p
    names=set(combined)|{p.name for p in root.glob('sample_*') if p.is_dir()}
    result=[]
    for sample in sorted(names):
        folder=root/sample
        candidates=[('Unicycler',folder/'unicycler/assembly.fasta','unicycler.done'),
                    ('Flye',folder/'flye/assembly.fasta','flye.done')]
        present=[entry for entry in candidates if entry[1].is_file() and entry[1].stat().st_size]
        # The combined FASTA is the one used in ANI; keep its identity as primary.
        source=combined.get(sample) or (present[0][1] if present else None)
        assembler=present[0][0] if len(present)==1 else 'Multiple outputs' if present else 'Not recorded'
        done=any((folder/entry[2]).is_file() for entry in present)
        row=dict(sample=sample,assembler=assembler,assembly_status='missing',in_combined=sample in combined,
                 stage_marker_present=done,source_file=str(source.relative_to(root)) if source else '',quast={})
        if source:
            try:
                row.update(fasta_stats(source));row['assembly_status']='completed' if done else 'FASTA present'
            except ValueError as exc:
                row['assembly_status']='invalid FASTA';row['assembly_error']=str(exc)
        row['quast']=quast_metrics(folder/'quast/report.tsv',sample)
        row['_path']=source
        result.append(row)
    return result


class NewickParser:
    def __init__(self,text): self.text=text;self.i=0
    def skip(self):
        while self.i<len(self.text):
            if self.text[self.i].isspace(): self.i+=1
            elif self.text[self.i]=='[':
                end=self.text.find(']',self.i)
                if end<0: raise ValueError('Unclosed Newick comment')
                self.i=end+1
            else: break
    def token(self):
        self.skip()
        if self.i<len(self.text) and self.text[self.i]=="'":
            self.i+=1;parts=[]
            while self.i<len(self.text):
                c=self.text[self.i];self.i+=1
                if c=="'":
                    if self.i<len(self.text) and self.text[self.i]=="'": parts.append("'");self.i+=1
                    else: return ''.join(parts)
                else: parts.append(c)
            raise ValueError('Unclosed quoted Newick label')
        start=self.i
        while self.i<len(self.text) and self.text[self.i] not in ':,();[': self.i+=1
        return self.text[start:self.i].strip()
    def node(self):
        self.skip();children=[]
        if self.i<len(self.text) and self.text[self.i]=='(':
            self.i+=1
            while True:
                children.append(self.node());self.skip()
                if self.i>=len(self.text): raise ValueError('Unclosed Newick node')
                c=self.text[self.i];self.i+=1
                if c==')': break
                if c!=',': raise ValueError('Invalid Newick separator')
        label=self.token();self.skip();length=None
        if self.i<len(self.text) and self.text[self.i]==':':
            self.i+=1;value=self.token();length=numeric(value)
            if length is None or length<0: raise ValueError('Invalid/negative Newick branch length: '+value)
        if not children and not label: raise ValueError('Empty Newick leaf')
        return dict(label=label,length=length,children=children)
    def parse(self):
        node=self.node();self.skip()
        if self.i>=len(self.text) or self.text[self.i]!=';': raise ValueError('Newick must end with semicolon')
        self.i+=1;self.skip()
        if self.i!=len(self.text): raise ValueError('Multiple Newick trees are not supported')
        return node


def leaves(node):
    if not node['children']: yield node
    for child in node['children']: yield from leaves(child)


def choose_tree(root,explicit=None):
    if explicit:
        path=explicit.resolve()
        if not path.is_file(): raise ValueError('Tree not found: '+str(path))
        return path
    paths=sorted({p.resolve() for p in root.glob('cohort_phylogen*/GToTree_output/*.tre')}|
                 {p.resolve() for p in root.glob('cohort_phylogen*/runs/attempt_*/*.tre')})
    valid=[]
    for path in paths:
        try:
            NewickParser(path.read_text()).parse();valid.append(path)
        except ValueError: pass
    if len(valid)>1: raise ValueError('Multiple trees found; choose one with --tree: '+', '.join(map(str,valid)))
    return valid[0] if valid else None


def tree_annotations(tree_path,tree,samples,nodes):
    mapping={}
    for sample in samples:
        for label in (sample['sample'],'sample_'+sample['sample']): mapping[norm(label)]=('sample',sample['sample'])
    phylo=next((p for p in tree_path.parents if (p/'selected_genomes.csv').is_file()),tree_path.parent.parent)
    for row in read_table(phylo/'selected_genomes.csv'):
        if row.get('tree_label') and row.get('id'):
            mapping[norm(row['tree_label'])]=(row.get('kind','reference'),row['id'])
        if row.get('input_fasta'):
            mapping[Path(row['input_fasta']).stem]=(row.get('kind','reference'),row.get('id',''))
    quality=read_table(tree_path.parent/'Genomes_summary_info.tsv','\t')
    metadata={}
    for row in quality:
        match=mapping.get(norm(row.get('label',''))) or mapping.get(norm(row.get('assembly_id','')))
        if match: metadata[match[1]]=row
    for leaf in leaves(tree):
        label=leaf['label'];found=mapping.get(norm(label))
        if found is None:
            match=re.search(r'GC[AF]_\d+\.\d+',label)
            found=('reference',match[0]) if match else ('unknown',label)
        kind,identifier=found
        row=nodes.get(identifier,{})
        leaf.update(kind=kind,id=identifier,organism=row.get('organism',''),strain=row.get('strain',''))
    return metadata


def read_ani(folder,samples):
    ids={s['sample'] for s in samples}
    sample_pairs=[];all_references={};seen=set()
    path=folder/'ani_pairs.csv'
    if not path.is_file(): return [],{},[]
    with path.open() as f:
        for row in csv.DictReader(f):
            a,b=row.get('genome_a'),row.get('genome_b')
            if not a or not b: raise ValueError('ani_pairs.csv lacks genome IDs')
            if a not in ids and b not in ids: continue
            key=tuple(sorted((a,b)))
            if key in seen: raise ValueError('Duplicate ANI pair: '+str(key))
            seen.add(key)
            ani=numeric(row.get('ani_used'));af=numeric(row.get('min_aligned_fraction'))
            accepted=row.get('status')=='accepted'
            if accepted and (ani is None or not 0<=ani<=100 or af is None or not 0<=af<=1):
                raise ValueError('Invalid accepted ANI record: '+str(key))
            item=dict(a=a,b=b,ani=ani if accepted else None,observed_ani=numeric(row.get('ani_mean')),
                      af=af,status=row.get('status','unknown'))
            if a in ids and b in ids: sample_pairs.append(item)
            elif accepted:
                sample,ref=(a,b) if a in ids else (b,a)
                if re.fullmatch(r'GC[AF]_\d+\.\d+',ref):
                    all_references.setdefault(sample,[]).append(dict(reference=ref,ani=ani,af=af))
    nearest={}
    for sample,values in all_references.items():
        values.sort(key=lambda r:(-r['ani'],-r['af'],r['reference']))
        nearest[sample]=values[0]
    return sample_pairs,nearest,sorted({r['reference'] for values in all_references.values() for r in values[:5]})


def similarity_groups(samples,pairs,threshold):
    """Complete linkage: every cross pair must be measured/accepted and exceed threshold."""
    names=sorted(s['sample'] if isinstance(s,dict) else s for s in samples)
    lookup={tuple(sorted((p['a'],p['b']))):p['ani'] for p in pairs if p['ani'] is not None}
    observed={x for pair in lookup for x in pair}
    groups=[[n] for n in names if n in observed]
    while True:
        choices=[]
        for i in range(len(groups)):
            for j in range(i+1,len(groups)):
                values=[lookup.get(tuple(sorted((a,b)))) for a in groups[i] for b in groups[j]]
                if all(v is not None and v>=threshold for v in values): choices.append((-min(values),tuple(groups[i]+groups[j]),i,j))
        if not choices: break
        _,_,i,j=min(choices);groups[i]=sorted(groups[i]+groups[j]);groups.pop(j)
    groups.sort(key=lambda g:(-len(g),g))
    return {name:'G%s' % (i+1) for i,g in enumerate(groups) for name in g}|{name:'Unresolved' for name in names if name not in observed}


def csv_bytes(rows,fields):
    output=io.StringIO(newline='');writer=csv.DictWriter(output,fieldnames=fields,extrasaction='ignore')
    writer.writeheader();writer.writerows(rows)
    return output.getvalue().encode()


def build(args):
    root=args.root.expanduser().resolve();out=args.outdir.expanduser().resolve()
    if not root.is_dir(): raise ValueError('Results root not found')
    if out==root or out in root.parents: raise ValueError('Use a separate dashboard output folder')
    owner=out/'bundle_manifest.json'
    if out.exists() and any(out.iterdir()) and not owner.is_file(): raise ValueError('Output folder is not an existing dashboard bundle; choose a new folder')
    samples=sample_files(root)
    if not samples: raise ValueError('No sample folders or combined assemblies found')
    pairs,nearest,refs=read_ani(root/'cohort_ani',samples)
    nodes={r['id']:r for r in read_table(root/'cohort_ani/network_nodes.csv')}
    tree_path=choose_tree(root,args.tree)
    tree=NewickParser(tree_path.read_text()).parse() if tree_path else None
    metadata=tree_annotations(tree_path,tree,samples,nodes) if tree else {}
    in_tree={n['id'] for n in leaves(tree) if n.get('kind')=='sample'} if tree else set()
    groups=similarity_groups(samples,pairs,args.group_ani)
    safe_names=[re.sub(r'[^A-Za-z0-9_.-]','_',s['sample']) for s in samples]
    if len(safe_names)!=len(set(safe_names)): raise ValueError('Sanitized sample filenames collide')
    files=[];out.mkdir(parents=True,exist_ok=True)
    for sample in samples:
        identifier=sample['sample'];source=sample.pop('_path')
        sample.update(ani_group=groups[identifier],in_tree=identifier in in_tree,
                      nearest_reference=nearest.get(identifier,{}).get('reference'),
                      nearest_ani=nearest.get(identifier,{}).get('ani'),
                      nearest_af=nearest.get(identifier,{}).get('af'))
        q=metadata.get(identifier,{})
        sample['marker_recovery_percent']=numeric(q.get('perc_comp'))
        sample['marker_redundancy_percent']=numeric(q.get('perc_redund'))
        if source and sample['assembly_status']!='invalid FASTA':
            safe=re.sub(r'[^A-Za-z0-9_.-]','_',identifier)
            target=out/'assemblies'/(safe+'.fa.gz');target.parent.mkdir(exist_ok=True)
            opener=gzip.open if str(source).endswith('.gz') else open
            with opener(source,'rb') as src,target.open('wb') as raw,gzip.GzipFile(fileobj=raw,mode='wb',mtime=0) as dest:
                shutil.copyfileobj(src,dest)
            sample['download']='assemblies/'+target.name;files.append(target)
    # Ensure filename sanitization cannot silently overwrite another sample.
    downloads=[s['download'] for s in samples if s.get('download')]
    if len(downloads)!=len(set(downloads)): raise ValueError('Sanitized sample filenames collide')
    fields=['sample','assembly_status','assembler','in_combined','in_tree','length_bp','contigs','n50','largest_contig','gc_percent',
            'circular_contigs','marker_recovery_percent','marker_redundancy_percent','ani_group','nearest_reference','nearest_ani','nearest_af']
    export=out/'data';export.mkdir(exist_ok=True)
    outputs={export/'assemblies.csv':csv_bytes(samples,fields),export/'sample_ani.csv':csv_bytes(pairs,['a','b','ani','observed_ani','af','status'])}
    if tree_path: outputs[export/'phylogeny.tre']=tree_path.read_bytes()
    for path,content in outputs.items(): path.write_bytes(content);files.append(path)
    accepted=[p['ani'] for p in pairs if p['ani'] is not None]
    payload=dict(version=VERSION,title=args.title,samples=samples,pairs=pairs,tree=tree,tree_file=tree_path.name if tree_path else None,
                 references=[dict(id=r,**{k:nodes.get(r,{}).get(k,'') for k in ['organism','strain','isolate']}) for r in refs],
                 threshold=args.group_ani,ani_min=min(accepted) if accepted else None,ani_max=max(accepted) if accepted else None,
                 accepted_pair_count=len(accepted),expected_pair_count=len(samples)*(len(samples)-1)//2,
                 generated=datetime.now(timezone.utc).isoformat())
    template=HTML_TEMPLATE.replace('__PAYLOAD__',json.dumps(payload,ensure_ascii=True).replace('<','\\u003c'))
    (out/'dashboard.html').write_text(template);files.append(out/'dashboard.html')
    manifest=dict(tool='spraynpray_dashboard',version=VERSION,files=[str(f.relative_to(out)) for f in files],
                  source_tree=str(tree_path) if tree_path else None,group_ani=args.group_ani)
    owner.write_text(json.dumps(manifest,indent=2)+'\n');files.append(owner)
    bundle=out.with_name(out.name+'.zip')
    with zipfile.ZipFile(bundle,'w',zipfile.ZIP_DEFLATED) as z:
        for path in files: z.write(path,out.name+'/'+str(path.relative_to(out)))
    print('Dashboard:',out/'dashboard.html');print('Portable ZIP:',bundle)
    print('Assembly FASTAs recovered: %s / %s; accepted sample ANI pairs: %s / %s; samples in tree: %s' % (
          sum(s['assembly_status'] in ['completed','FASTA present'] for s in samples),len(samples),len(accepted),payload['expected_pair_count'],len(in_tree)))
    return payload


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True,help='hybrid_results folder')
    p.add_argument('--outdir',type=Path,required=True,help='Portable dashboard output folder')
    p.add_argument('--tree',type=Path,help='Explicit GToTree Newick file; auto-detects when exactly one exists')
    p.add_argument('--title',default='Isolate genome comparison')
    p.add_argument('--group-ani',type=float,default=99,help='Exploratory complete-linkage similarity threshold, percent; not a lineage/transmission cutoff')
    args=p.parse_args(argv)
    if not math.isfinite(args.group_ani) or not 0<=args.group_ani<=100: p.error('--group-ani must be 0–100')
    try: build(args)
    except (OSError,ValueError,KeyError,csv.Error) as exc: p.exit(1,str(exc)+'\n')
    return 0

HTML_TEMPLATE = r'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Isolate genome comparison</title>
<style>
:root{font-family:system-ui,-apple-system,sans-serif;color:#192c3d;background:#f3f6f9;font-size:15px;--accent:#087e83;--muted:#5b6d7d;--line:#dae2e9}*{box-sizing:border-box}body{margin:0}header{background:#142e42;color:white;padding:28px max(24px,calc((100vw - 1360px)/2))}header small{letter-spacing:.12em;text-transform:uppercase;color:#9cdbd7}h1{font-size:30px;line-height:1.2;margin:9px 0}h2{font-size:21px;margin:0 0 15px}h3{font-size:16px}p{line-height:1.6}header p{margin:8px 0;color:#cbdae5}main{max-width:1408px;margin:auto;padding:24px}nav{display:flex;gap:8px;flex-wrap:wrap;margin-bottom:24px}button,.download{border:1px solid #bacbd6;border-radius:6px;padding:9px 13px;background:white;color:#173e51;cursor:pointer;font:inherit;text-decoration:none;display:inline-block}button:hover,.download:hover{background:#e9f4f3}button.active{background:var(--accent);border-color:var(--accent);color:white}input,select{font:inherit;border:1px solid #b5c4cf;border-radius:5px;padding:8px;background:white;color:inherit}input[type=checkbox]{accent-color:var(--accent)}.panel{display:none}.panel.active{display:block}.card{background:white;border:1px solid var(--line);border-radius:10px;padding:24px;margin-bottom:20px}.cards{display:grid;grid-template-columns:repeat(4,1fr);gap:16px}.metric b{font-size:29px;display:block;margin:7px 0}.muted,.note{color:var(--muted);font-size:13px}.note{line-height:1.65}.banner{border-left:5px solid var(--accent)}.toolbar{display:flex;gap:12px;align-items:center;flex-wrap:wrap;margin:14px 0}.scroll{overflow:auto;max-height:640px}table{border-collapse:collapse;white-space:nowrap;width:100%;font-size:13px}th,td{text-align:left;padding:11px 12px;border-bottom:1px solid #e3e9ee}th{background:#eef3f7;position:sticky;top:0;z-index:1}th button{font-size:12px;padding:0;border:0;background:transparent!important}tr.selectable{cursor:pointer}tr.selectable:hover{background:#edf8f7}.pill{display:inline-block;padding:3px 7px;border-radius:4px;background:#e3f1ed;color:#12614d}.warn{background:#fff0d6;color:#7c5010}.legend{display:flex;gap:15px;flex-wrap:wrap;font-size:13px;align-items:center}.swatch{display:inline-block;width:11px;height:11px;border-radius:3px;margin-right:5px;vertical-align:middle}.detail{background:#f7fafc;border-radius:6px;padding:17px;margin-top:15px}.detail dl{display:grid;grid-template-columns:max-content 1fr;gap:6px 20px;font-size:13px}.detail dd{margin:0;overflow-wrap:anywhere}#treeViewport{overflow:auto;max-height:740px;border:1px solid var(--line);border-radius:6px;background:white}#heatmap{display:block;cursor:crosshair}#tooltip{position:fixed;display:none;background:#152e42;color:white;border-radius:5px;padding:9px 12px;pointer-events:none;font-size:13px;z-index:10;max-width:420px}a{color:#06777f}footer{font-size:12px;color:var(--muted);padding:0 24px 25px;max-width:1408px;margin:auto}svg text{font-family:system-ui,sans-serif}li{line-height:1.7}button:focus-visible,a:focus-visible,input:focus-visible{outline:3px solid #e6a340;outline-offset:3px}@media(max-width:800px){.cards{grid-template-columns:repeat(2,1fr)}main{padding:15px}.card{padding:17px}h1{font-size:25px}}
</style></head><body>
<header><small>Middle Author Bioinformatics · cohort review</small><h1 id="title"></h1><p>Assembly outcomes, genome similarity and phylogenomic context</p></header>
<main>
<nav aria-label="Dashboard sections"><button class="active" data-tab="overview">Overview</button><button data-tab="assemblies">Assemblies</button><button data-tab="ani">ANI comparisons</button><button data-tab="phylogeny">Phylogeny</button></nav>
<section id="overview" class="panel active">
<div class="cards"><div class="card metric">Assemblies recovered<b id="assemblyCount"></b><span class="muted">Valid final FASTAs / sample IDs found</span></div><div class="card metric">Sample ANI coverage<b id="pairCount"></b><span class="muted">Accepted / possible sample pairs</span></div><div class="card metric">Samples in phylogeny<b id="treeCount"></b><span class="muted">Matched to native tree labels</span></div><div class="card metric">Similarity groups<b id="groupCount"></b><span class="muted" id="groupCaption"></span></div></div>
<div class="card banner"><h2 id="finding"></h2><p id="findingText"></p><p class="note">These are exploratory similarity groups, not validated strain or transmission assignments. ANI and a marker-gene tree cannot by themselves establish a single evolving lineage, isolate provenance or direction of transmission.</p></div>
<div class="card"><h2>Follow-up: reconcile isolate origins</h2><p>Review the ANI groups and phylogeny alongside the collection history before treating the cohort as one evolving lineage.</p><ul><li>Match sample IDs to isolate source, collection date, location and culture history.</li><li>Review samples outside the main groups and assemblies absent from ANI or the tree.</li><li>For closely related isolates, decide whether read-level QC, reference-based SNP analysis and recombination assessment are needed.</li></ul></div>
<div class="card"><h2>Analysis coverage and quality</h2><div id="coverageNotes"></div><p class="note">Lengths, contig counts, N50 and GC are computed from the packaged final FASTAs. QUAST report values are available in each sample’s details. GToTree marker recovery and redundancy describe the selected marker set; they are not CheckM completeness or contamination estimates. No CheckM assessment is inferred.</p><div class="toolbar"><a class="download" href="data/assemblies.csv" download>Assembly summary CSV</a><a class="download" href="data/sample_ani.csv" download>Sample ANI CSV</a><button id="exportGroups">Current groups CSV</button></div></div>
</section>
<section id="assemblies" class="panel"><div class="card"><h2>Assemblies and analysis status</h2><p class="note">“Completed” means a valid final FASTA and an assembler completion marker were found. It is not a biological quality pass. Click a row for QUAST, marker metrics, nearest reference and an assembly download.</p><div class="toolbar"><label>Search <input id="sampleSearch" type="search" placeholder="Sample ID"></label><label>Show <select id="statusFilter"><option value="all">All samples</option><option value="outside">Absent from combined analysis</option><option value="missing">Missing / invalid assemblies</option><option value="notree">Absent from phylogeny</option></select></label><button id="clearFilters">Clear filters</button><span id="rowCount" class="muted"></span></div><div class="scroll"><table><thead><tr id="assemblyHead"></tr></thead><tbody id="assemblyBody"></tbody></table></div><div class="detail" id="sampleDetail" aria-live="polite">Select a sample to view its details.</div></div></section>
<section id="ani" class="panel"><div class="card"><h2>Genome-wide similarity among input samples</h2><div class="toolbar"><label>Grouping threshold <input id="threshold" type="number" min="0" max="100" step="0.01" style="width:90px"> % ANI</label><button id="applyThreshold">Update groups</button><span id="thresholdError" role="alert"></span></div><p class="note">Complete linkage: groups merge only when every cross-pair has an accepted ANI at or above the threshold. Unmeasured or rejected pairs prevent a merge; samples with no accepted sample comparison are “Unresolved.” The default 99% is an exploration setting, not a biological lineage cutoff. Only sample–sample comparisons determine these groups.</p><div class="legend"><span><i class="swatch" style="background:#0b6170"></i>Higher ANI</span><span><i class="swatch" style="background:#e7f3ed"></i>Lower ANI</span><span><i class="swatch" style="background:#dfe3e8"></i>Unknown / rejected</span><span id="aniRange"></span></div><div class="scroll" style="max-height:900px;margin-top:18px"><canvas id="heatmap" aria-label="Pairwise sample ANI heatmap"></canvas></div><div id="pairDetail" class="detail" aria-live="polite">Hover or click a cell for exact ANI and aligned fraction.</div><div id="groupLegend" class="legend" style="margin-top:18px"></div></div></section>
<section id="phylogeny" class="panel"><div class="card"><h2>Marker-gene phylogeny</h2><p class="note">Input samples are shown by default. Including references displays all leaves in the supplied tree. Labels have dedicated rows. Branch lengths follow the native tree; displayed rooting does not establish an ancestor or temporal direction.</p><div class="toolbar"><label><input id="includeReferences" type="checkbox"> Include reference genomes</label><label>Find label <input id="treeSearch" type="search"></label><button id="zoomIn">Zoom in</button><button id="zoomOut">Zoom out</button><button id="resetTree">Reset</button><button id="saveTree">Download displayed SVG</button><a class="download" id="newickLink" href="data/phylogeny.tre" download>Native Newick</a></div><p id="treeNote" class="note"></p><div id="treeViewport"></div><div id="treeDetail" class="detail" aria-live="polite">Select a leaf for genome information.</div></div></section>
</main><footer id="footer"></footer><div id="tooltip" role="tooltip"></div>
<script id="payload" type="application/json">__PAYLOAD__</script>
<script>
'use strict';
const D=JSON.parse(document.getElementById('payload').textContent);
const $=id=>document.getElementById(id), esc=x=>String(x??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const fmt=(v,d=0)=>v===null||v===undefined?'Not available':Number(v).toLocaleString(undefined,{maximumFractionDigits:d,minimumFractionDigits:d});
const good=s=>['completed','FASTA present'].includes(s.assembly_status), key=(a,b)=>JSON.stringify([a,b].sort());
const lookup=new Map(D.pairs.map(p=>[key(p.a,p.b),p])), samples=new Map(D.samples.map(s=>[s.sample,s]));
const colors=['#087f8c','#8b5fbf','#d17c22','#438655','#b34e68','#3d70b5','#93784a','#6c7e8e'];
let groups={},threshold=D.threshold,sortField='sample',sortUp=true,treeZoom=1;
function color(group){return group==='Unresolved'?'#8996a3':colors[(Number(group.slice(1))-1)%colors.length];}
function leafList(node){return node.children.length?node.children.flatMap(leafList):[node];}
const treeLeaves=D.tree?leafList(D.tree):[];
const order=[...new Set([...treeLeaves.filter(l=>l.kind==='sample').map(l=>l.id),...D.samples.map(s=>s.sample)])].filter(s=>samples.has(s));
function makeGroups(){
 const accepted=D.pairs.filter(p=>p.ani!==null),observed=new Set(accepted.flatMap(p=>[p.a,p.b]));
 let clusters=[...samples.keys()].sort().filter(n=>observed.has(n)).map(n=>[n]);
 const cmp=(a,b)=>a<b?-1:a>b?1:0,seq=(a,b)=>{for(let i=0;i<Math.min(a.length,b.length);i++){const c=cmp(a[i],b[i]);if(c)return c;}return a.length-b.length;};
 while(true){let best=null;for(let i=0;i<clusters.length;i++)for(let j=i+1;j<clusters.length;j++){
  const vals=clusters[i].flatMap(a=>clusters[j].map(b=>lookup.get(key(a,b))?.ani??null));
  if(vals.every(v=>v!==null&&v>=threshold)){const proposal={i,j,min:Math.min(...vals),tie:[...clusters[i],...clusters[j]]};if(!best||proposal.min>best.min||(proposal.min===best.min&&seq(proposal.tie,best.tie)<0))best=proposal;}
 }if(!best)break;clusters[best.i]=[...clusters[best.i],...clusters[best.j]].sort();clusters.splice(best.j,1);}
 clusters.sort((a,b)=>b.length-a.length||seq(a,b));groups={};clusters.forEach((g,i)=>g.forEach(s=>groups[s]='G'+(i+1)));for(const s of samples.keys())if(!groups[s])groups[s]='Unresolved';
 return clusters;
}
function summary(){
 const clusters=makeGroups(),accepted=D.pairs.filter(p=>p.ani!==null),below=accepted.filter(p=>p.ani<threshold),unknown=D.expected_pair_count-accepted.length;
 $('assemblyCount').textContent=D.samples.filter(good).length+' / '+D.samples.length;
 $('pairCount').textContent=accepted.length+' / '+D.expected_pair_count;
 $('treeCount').textContent=D.samples.filter(s=>s.in_tree).length+' / '+D.samples.length;
 $('groupCount').textContent=clusters.length;$('groupCaption').textContent='At '+threshold+'% ANI; '+Object.values(groups).filter(g=>g==='Unresolved').length+' unresolved';
 if(below.length){$('finding').textContent='Evidence of distinct genomic similarity groups';$('findingText').textContent=below.length+' measured sample pair'+(below.length===1?' is':'s are')+' below '+threshold+'% ANI. Accepted sample ANI spans '+fmt(D.ani_min,3)+'–'+fmt(D.ani_max,3)+'%. This pattern is consistent with multiple genomic backgrounds at this threshold and warrants reviewing isolate origins. '+(unknown?unknown+' possible pair comparisons remain unknown or rejected.':'All possible sample pairs have accepted ANI.');}
 else {$('finding').textContent=accepted.length?'Available ANI does not establish separate lineages':'Sample ANI is not yet available';$('findingText').textContent=accepted.length?'No accepted sample pair is below the selected threshold. Separate groups can still arise from missing comparisons; group counts alone are not evidence of distinct lineages.':'No accepted sample–sample comparisons were found. Assembly availability and tree membership are shown, but a lineage conclusion cannot be drawn from missing ANI.';}
 const outside=D.samples.filter(s=>!s.in_combined).map(s=>s.sample),missing=D.samples.filter(s=>!good(s)).map(s=>s.sample),notree=D.samples.filter(s=>!s.in_tree).map(s=>s.sample);
 $('coverageNotes').innerHTML=[['Absent from combined assembly analysis',outside],['Missing or invalid final FASTA',missing],['Not present among matched phylogeny leaves',notree]].map(([label,names])=>'<p><strong>'+label+':</strong> '+(names.length?names.map(esc).join(', '):'None')+'</p>').join('');
 $('groupLegend').innerHTML=clusters.map((c,i)=>'<span><i class="swatch" style="background:'+color('G'+(i+1))+'"></i>G'+(i+1)+': '+c.length+' sample'+(c.length===1?'':'s')+'</span>').join('');
 $('aniRange').textContent=accepted.length?'Color range: '+fmt(D.ani_min,3)+'–100% ANI':'No accepted pairs';
 renderAssemblies();renderHeatmap();renderTree();
}
const columns=[['sample','Sample'],['assembly_status','Assembly'],['assembler','Assembler'],['length_bp','Length (bp)'],['contigs','Contigs'],['n50','N50 (bp)'],['gc_percent','GC %'],['ani_group','Group'],['in_combined','Combined'],['in_tree','Tree'],['nearest_ani','Nearest ref ANI %']];
function renderAssemblies(){
 $('assemblyHead').innerHTML=columns.map(([id,label])=>'<th><button data-sort="'+id+'">'+label+(id===sortField?(sortUp?' ↑':' ↓'):'')+'</button></th>').join('');
 let rows=D.samples.filter(s=>s.sample.toLowerCase().includes($('sampleSearch').value.toLowerCase()));const mode=$('statusFilter').value;
 rows=rows.filter(s=>mode==='all'||mode==='outside'&&!s.in_combined||mode==='missing'&&!good(s)||mode==='notree'&&!s.in_tree);
 const val=s=>sortField==='ani_group'?groups[s.sample]:s[sortField];rows.sort((a,b)=>{let x=val(a),y=val(b);if(x==null)return 1;if(y==null)return -1;return (typeof x==='string'?x.localeCompare(y):x-y)*(sortUp?1:-1);});
 $('assemblyBody').innerHTML=rows.map(s=>'<tr class="selectable" tabindex="0" data-sample="'+esc(s.sample)+'">'+columns.map(([id])=>'<td>'+(id==='assembly_status'?'<span class="pill '+(good(s)?'':'warn')+'">'+esc(s[id])+'</span>':id==='ani_group'?'<span style="color:'+color(groups[s.sample])+'">'+groups[s.sample]+'</span>':typeof s[id]==='boolean'?(s[id]?'Yes':'No'):typeof s[id]==='number'?fmt(s[id],id==='gc_percent'?2:id==='nearest_ani'?4:0):esc(s[id]??'—'))+'</td>').join('')+'</tr>').join('');
 $('rowCount').textContent=rows.length+' samples shown';
}
function details(s){
 const metrics=[['Assembly',s.assembly_status],['Assembler',s.assembler],['Combined assembly set',s.in_combined?'Included':'Not included'],['Phylogeny',s.in_tree?'Present':'Not present'],['Similarity group',groups[s.sample]],['Length (bp)',fmt(s.length_bp)],['Contigs',fmt(s.contigs)],['N50 (bp)',fmt(s.n50)],['Largest contig (bp)',fmt(s.largest_contig)],['GC (%)',fmt(s.gc_percent,2)],['Ambiguous bases',fmt(s.ambiguous_bp)],['Circular contigs (FASTA annotations)',fmt(s.circular_contigs)],['GToTree marker recovery (%)',fmt(s.marker_recovery_percent,2)],['GToTree marker redundancy (%)',fmt(s.marker_redundancy_percent,2)],['Nearest accepted reference',s.nearest_reference??'Not available'],['Reciprocal genome ANI (%)',fmt(s.nearest_ani,4)],['Minimum reciprocal aligned fraction',fmt(s.nearest_af,3)]];
 const ref=D.references.find(r=>r.id===s.nearest_reference);if(ref)for(const k of ['organism','strain','isolate'])metrics.push(['Reference '+k,ref[k]||'Not provided']);
 return '<h3>'+esc(s.sample)+'</h3><dl>'+metrics.map(([k,v])=>'<dt>'+esc(k)+'</dt><dd>'+esc(v)+'</dd>').join('')+'</dl>'+(s.download?'<a class="download" href="'+esc(s.download)+'" download>Download assembly FASTA (.gz)</a>':'')+(s.assembly_error?'<p>'+esc(s.assembly_error)+'</p>':'')+'<h3>QUAST report</h3>'+(Object.keys(s.quast).length?'<dl>'+Object.entries(s.quast).map(([k,v])=>'<dt>'+esc(k)+'</dt><dd>'+esc(v)+'</dd>').join('')+'</dl>':'<p class="note">No unambiguous per-sample QUAST report was found.</p>');
}
function selectSample(id,target='sampleDetail'){const s=samples.get(id);if(s)$(target).innerHTML=details(s);}
const cell=18,margin=180;let heatN=0;
function heatColor(ani){const lower=Math.min(D.ani_min??95,99.99),t=Math.max(0,Math.min(1,(ani-lower)/(100-lower)));return 'rgb('+[231+(11-231)*t,243+(97-243)*t,237+(112-237)*t].map(Math.round).join(',')+')';}
function renderHeatmap(){
 const canvas=$('heatmap'),ctx=canvas.getContext('2d');heatN=order.length;const size=margin+cell*heatN+20,dpr=window.devicePixelRatio||1;canvas.width=size*dpr;canvas.height=size*dpr;canvas.style.width=size+'px';canvas.style.height=size+'px';ctx.scale(dpr,dpr);ctx.fillStyle='#fff';ctx.fillRect(0,0,size,size);ctx.font='11px system-ui';
 order.forEach((a,i)=>{ctx.fillStyle=color(groups[a]);ctx.fillRect(margin-8,margin+i*cell,5,cell-1);ctx.fillStyle='#203d51';ctx.textAlign='right';ctx.fillText(a,margin-14,margin+i*cell+13);ctx.save();ctx.translate(margin+i*cell+13,margin-13);ctx.rotate(-Math.PI/2);ctx.textAlign='left';ctx.fillText(a,0,0);ctx.restore();
 order.forEach((b,j)=>{const ani=i===j?(good(samples.get(a))?100:null):lookup.get(key(a,b))?.ani??null;ctx.fillStyle=ani===null?'#dfe3e8':heatColor(ani);ctx.fillRect(margin+j*cell,margin+i*cell,cell-1,cell-1);});});
}
function pairText(a,b){if(a===b)return a+' · self comparison'+(good(samples.get(a))?' (100% by definition)':' · no valid assembly');const p=lookup.get(key(a,b));return a+' ↔ '+b+' · '+(p?.ani!=null?fmt(p.ani,5)+'% ANI · min aligned fraction '+fmt(p.af,3):'ANI unknown / rejected ('+(p?.status??'not measured')+')');}
function heatEvent(e,click){const r=$('heatmap').getBoundingClientRect(),j=Math.floor((e.clientX-r.left-margin)/cell),i=Math.floor((e.clientY-r.top-margin)/cell);if(i<0||j<0||i>=heatN||j>=heatN){$('tooltip').style.display='none';return;}const message=pairText(order[i],order[j]);if(click)$('pairDetail').textContent=message;else{const tip=$('tooltip');tip.textContent=message;tip.style.display='block';tip.style.left=Math.max(8,Math.min(innerWidth-420,e.clientX+15))+'px';tip.style.top=Math.min(innerHeight-80,e.clientY+15)+'px';}}
function prune(node,include){if(!node.children.length)return include||node.kind==='sample'?{...node}:null;const children=node.children.map(c=>prune(c,include)).filter(Boolean);if(!children.length)return null;if(children.length===1){const c=children[0];return {...c,length:node.length===null||c.length===null?null:node.length+c.length};}return {...node,children};}
function renderTree(){
 const viewport=$('treeViewport');viewport.replaceChildren();if(!D.tree){$('treeNote').textContent='No phylogeny file found. Provide --tree when generating this dashboard.';$('newickLink').hidden=true;$('saveTree').disabled=true;return;}
 const tree=prune(D.tree,$('includeReferences').checked);if(!tree){$('treeNote').textContent='No sample labels could be matched. Include references to inspect all native labels.';return;}
 const ls=leafList(tree),allNodes=[];function walk(n,d){allNodes.push(n);n.x=d;n.children.forEach(c=>walk(c,d+(c.length??1)));}walk(tree,0);
 const scaled=allNodes.filter(n=>n!==tree).every(n=>n.length!==null);if(!scaled){function topo(n,d){n.x=d;n.children.forEach(c=>topo(c,d+1));}topo(tree,0);}
 const row=23*treeZoom,width=1100*treeZoom,height=Math.max(140,(ls.length+2)*row),left=25,max=Math.max(...allNodes.map(n=>n.x),1e-12),span=650*treeZoom;
 ls.forEach((n,i)=>n.y=(i+1)*row);function position(n){if(n.children.length){n.children.forEach(position);n.y=n.children.reduce((s,c)=>s+c.y,0)/n.children.length;}}position(tree);
 const NS='http://www.w3.org/2000/svg',svg=document.createElementNS(NS,'svg');svg.setAttribute('xmlns',NS);svg.setAttribute('width',width);svg.setAttribute('height',height);svg.setAttribute('viewBox','0 0 '+width+' '+height);svg.setAttribute('role','img');svg.setAttribute('aria-label','Phylogenomic tree');
 const element=(tag,attrs,text)=>{const el=document.createElementNS(NS,tag);for(const [k,v] of Object.entries(attrs))el.setAttribute(k,v);if(text!==undefined)el.textContent=text;svg.appendChild(el);return el;},x=n=>left+n.x/max*span;
 function draw(n){if(n.children.length){element('line',{x1:x(n),x2:x(n),y1:Math.min(...n.children.map(c=>c.y)),y2:Math.max(...n.children.map(c=>c.y)),stroke:'#7b8d9b','stroke-width':1});for(const c of n.children){const branch=element('line',{x1:x(n),x2:x(c),y1:c.y,y2:c.y,stroke:'#7b8d9b','stroke-width':1});const title=document.createElementNS(NS,'title');title.textContent='Branch length: '+(c.length??'not recorded')+(c.children.length&&c.label?' · internal label/support: '+c.label:'');branch.appendChild(title);draw(c);}}
 else{const query=$('treeSearch').value.toLowerCase(),matched=query&&(n.id+' '+n.label+' '+n.organism+' '+n.strain).toLowerCase().includes(query),g=n.kind==='sample'?groups[n.id]:null;const label=n.kind==='sample'?n.id:n.id+(n.organism?' · '+n.organism:'')+(n.strain?' · '+n.strain:'');const circle=element('circle',{cx:x(n),cy:n.y,r:matched?5:3.5,fill:n.kind==='sample'?'#111':'#aaa'});const text=element('text',{x:x(n)+9,y:n.y+4,'font-size':12*treeZoom,fill:matched?'#ad5e00':'#233d50','font-weight':matched?'bold':'normal'},label);if(g)element('rect',{x:x(n)-5,y:n.y-8,width:3,height:16,fill:color(g)});
 [circle,text].forEach(el=>{el.style.cursor='pointer';el.addEventListener('click',()=>{if(n.kind==='sample')selectSample(n.id,'treeDetail');else $('treeDetail').innerHTML='<h3>'+esc(n.id)+'</h3><p>'+esc(n.organism||'Organism not recorded')+'</p><p>Strain: '+esc(n.strain||'Not provided')+'</p><p class="note">Native label: '+esc(n.label)+'</p>';});});}}
 draw(tree);const scale=scaled?max/5:1,scalePx=scale/max*span;const yy=height-10;element('line',{x1:left,x2:left+scalePx,y1:yy,y2:yy,stroke:'#243e50'});element('text',{x:left+scalePx+8,y:yy+4,'font-size':11},scaled?scale.toPrecision(3)+' substitutions/site':'1 branch');viewport.appendChild(svg);
 // Long reference labels get real horizontal space rather than clipping or overlapping.
 try{const box=svg.getBBox();svg.setAttribute('width',Math.max(width,box.x+box.width+24));svg.setAttribute('viewBox','0 0 '+Math.max(width,box.x+box.width+24)+' '+height);}catch(e){}
 $('treeNote').textContent=ls.length+' displayed leaves · '+(scaled?'Native branch lengths; collapsed paths retain their summed length.':'Some branch lengths are missing; this view uses topological spacing.')+' Native file: '+D.tree_file+'.';
}
function download(name,text,type){const a=document.createElement('a'),url=URL.createObjectURL(new Blob([text],{type}));a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);}
document.querySelectorAll('[data-tab]').forEach(b=>b.addEventListener('click',()=>{document.querySelectorAll('[data-tab]').forEach(x=>x.classList.toggle('active',x===b));document.querySelectorAll('.panel').forEach(s=>s.classList.toggle('active',s.id===b.dataset.tab));if(b.dataset.tab==='phylogeny')renderTree();}));
$('assemblyHead').addEventListener('click',e=>{const b=e.target.closest('[data-sort]');if(!b)return;if(sortField===b.dataset.sort)sortUp=!sortUp;else{sortField=b.dataset.sort;sortUp=true;}renderAssemblies();});
$('assemblyBody').addEventListener('click',e=>{const row=e.target.closest('[data-sample]');if(row)selectSample(row.dataset.sample);});$('assemblyBody').addEventListener('keydown',e=>{if(e.key==='Enter'){const row=e.target.closest('[data-sample]');if(row)selectSample(row.dataset.sample);}});
$('sampleSearch').addEventListener('input',renderAssemblies);$('statusFilter').addEventListener('change',renderAssemblies);$('clearFilters').addEventListener('click',()=>{$('sampleSearch').value='';$('statusFilter').value='all';renderAssemblies();});
$('threshold').value=threshold;$('applyThreshold').addEventListener('click',()=>{const value=Number($('threshold').value);if($('threshold').value===''||!Number.isFinite(value)||value<0||value>100){$('thresholdError').textContent='Enter a threshold from 0 to 100.';return;}$('thresholdError').textContent='';threshold=value;summary();});
$('heatmap').addEventListener('mousemove',e=>heatEvent(e,false));$('heatmap').addEventListener('click',e=>heatEvent(e,true));$('heatmap').addEventListener('mouseleave',()=>$('tooltip').style.display='none');
$('includeReferences').addEventListener('change',renderTree);$('treeSearch').addEventListener('input',renderTree);$('zoomIn').addEventListener('click',()=>{treeZoom=Math.min(3,treeZoom*1.25);renderTree();});$('zoomOut').addEventListener('click',()=>{treeZoom=Math.max(.65,treeZoom/1.25);renderTree();});$('resetTree').addEventListener('click',()=>{treeZoom=1;$('treeSearch').value='';$('includeReferences').checked=false;renderTree();});
$('saveTree').addEventListener('click',()=>{const svg=$('treeViewport').querySelector('svg');if(svg)download('phylogeny.svg',new XMLSerializer().serializeToString(svg),'image/svg+xml');});
$('exportGroups').addEventListener('click',()=>{const quote=x=>'"'+String(x).replace(/"/g,'""')+'"';download('similarity_groups.csv','sample,group,threshold_ani_percent\n'+D.samples.map(s=>[s.sample,groups[s.sample],threshold].map(quote).join(',')).join('\n')+'\n','text/csv');});
$('title').textContent=D.title;document.title=D.title;$('footer').textContent='Generated '+new Date(D.generated).toLocaleString()+' · SprayNPray dashboard '+D.version+' · Offline bundle: keep dashboard.html, data/ and assemblies/ together. No analysis is rerun by opening this page.';summary();
</script></body></html>
'''

if __name__=='__main__':
    raise SystemExit(main())
