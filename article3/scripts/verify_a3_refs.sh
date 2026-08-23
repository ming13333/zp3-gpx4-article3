#!/bin/bash
# Crossref DOI 批量核验（curl 版）
OUT="C:/D/workbuddy/科研/细胞外GPX4免疫抑制/article3/reports/a3_crossref_raw.jsonl"
> "$OUT"
PY="C:/Users/lm962/.workbuddy/binaries/python/versions/3.13.12/python.exe"

parse_one() {
  local n="$1"; local doi="$2"
  local body
  body=$(curl -s --max-time 25 "https://api.crossref.org/works/${doi}")
  echo "${body}" | "$PY" -c "
import sys, json
n = int(sys.argv[1]); doi = sys.argv[2]
raw = sys.stdin.read()
try:
    j = json.loads(raw)['message']
except Exception:
    print(json.dumps({'n':n,'doi':doi,'status':'FAIL','result':raw[:150]}, ensure_ascii=False))
    sys.exit()
title = ' '.join(j.get('title') or [''])[:200]
ct = (j.get('container-title') or [''])
container = ct[0] if ct else ''
year = None
for f in ['published-print','published-online','issued','created']:
    d = j.get(f)
    if d and d.get('date-parts') and d['date-parts'][0][0]:
        year = d['date-parts'][0][0]; break
vol = j.get('volume',''); issue = j.get('issue','')
pg = j.get('page','') or j.get('article-number','')
auths = [(a.get('family',''), a.get('given','')) for a in (j.get('author') or [])][:6]
print(json.dumps({'n':n,'doi':doi,'status':'OK','year':year,'container':container,
 'volume':vol,'issue':issue,'page':pg,'title':title,'authors':auths}, ensure_ascii=False))
" "$n" "$doi" >> "$OUT"
}

while IFS=$'\t' read -r n doi; do
  [ -z "$n" ] && continue
  echo "[${n}] ${doi}"
  parse_one "$n" "$doi"
  sleep 0.4
done << 'REFS'
1	10.1002/bies.202400248
2	10.1093/bib/bbx122
3	10.1016/j.cell.2024.02.013
4	10.3390/ijms27135866
5	10.3390/biomedicines12122850
6	10.1002/ijc.35098
7	10.3389/fonc.2023.1233039
8	10.1016/0012-1606(80)90371-1
9	10.1016/j.mce.2021.111502
10	10.1016/j.tice.2018.11.001
11	10.1186/s12957-026-04245-2
12	10.2174/1386207325666221010112601
13	10.2174/0109298665350171241204153202
14	10.1038/s41598-021-86888-7
15	10.1016/j.cell.2025.12.002
16	10.1080/2162402X.2026.2624244
17	10.1038/s41392-020-00356-8
18	10.32604/or.2026.082613
19	10.1002/0471142905.hg1116s87
20	10.1093/nar/gkv1288
21	10.1038/nature07385
22	10.1038/ng.2653
23	10.1038/ng.3969
24	10.1073/pnas.0506580102
25	10.18637/jss.v067.i01
26	10.3389/fonc.2026.1888793
27	10.1016/j.prp.2026.156634
28	10.3389/fimmu.2026.1857659
29	10.1093/neuonc/noag151
30	10.1016/j.expneurol.2026.115888
31	10.1097/SCS.0000000000013133
32	10.1016/j.mtbio.2026.103491
33	10.1186/s12974-026-03964-3
34	10.1007/s11481-026-10302-0
35	10.1016/j.immuni.2026.06.024
36	10.1016/j.intimp.2026.117206
37	10.1093/neuonc/nou307
38	10.1186/s12935-026-04384-2
39	10.2174/0115680096474313260623050646
40	10.1016/j.compbiolchem.2026.109165
41	10.3389/fimmu.2026.1820802
42	10.2174/0115748928461929260701074511
43	10.1016/j.smim.2026.102049
44	10.1038/s41420-026-03214-8
45	10.3389/fimmu.2026.1848122
46	10.1007/s10142-026-01973-2
47	10.1186/s12885-026-16321-7
48	10.55782/z4j0rp75
49	10.1016/j.nbd.2026.107530
50	10.1016/j.medj.2026.101228
51	10.1186/s12967-026-08555-7
52	10.3390/genes17060610
53	10.1371/journal.pone.0351849
54	10.32604/or.2026.079221
55	10.1186/s42466-026-00515-2
56	10.1093/noajnl/vdag182
57	10.1155/ijog/2045937
58	10.1016/j.canlet.2026.218716
59	10.1016/j.compbiolchem.2026.109202
60	10.1016/j.intimp.2026.117056
61	10.3389/fimmu.2026.1825000
62	10.1007/s10142-026-01918-9
63	10.1016/j.brainresbull.2026.112053
64	10.3389/fimmu.2026.1871988
65	10.1016/j.canlet.2026.218718
66	10.1186/s12967-026-08497-0
67	10.1186/s13062-026-00889-y
68	10.3390/cancers18132092
69	10.1080/15384101.2026.2703934
70	10.1038/s41388-026-03844-3
71	10.1016/j.clim.2026.110755
72	10.1186/s12885-026-16332-4
73	10.1016/j.canlet.2026.218717
76	10.1172/jci.insight.85841
77	10.1038/s41556-018-0118-z
78	10.1038/nbt.3519
REFS
echo "DONE $(wc -l < "$OUT") records"
