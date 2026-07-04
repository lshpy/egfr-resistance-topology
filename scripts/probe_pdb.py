"""EGFR 내성 캐스케이드 후보 PDB 구조의 실제 상태를 프로그램으로 판별.
각 구조에서 결합 리간드(HETATM)와 키 잔기(790/797/858) 정체를 읽는다.
"""
import warnings, urllib.request, os
warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser

CANDIDATES = [
    # (추정 상태, PDB ID)
    ("WT/erlotinib", "1M17"), ("WT/lapatinib", "1XKK"), ("WT/gefitinib", "2ITY"),
    ("WT/gefitinib", "2ITZ"), ("WT/AEE788", "2ITX"), ("WT", "2J6M"),
    ("WT/TAK285", "3POZ"), ("WT", "4LQM"), ("WT/afatinib", "4G5J"),
    ("T790M", "3W2Q"), ("T790M", "3W2R"), ("T790M/AEE788", "3IKA"),
    ("T790M", "2JIU"), ("T790M/afatinib", "4G5P"), ("T790M", "4I22"),
    ("T790M", "4I23"), ("T790M/osimertinib?", "4ZAU"), ("T790M/osi", "6JX0"),
    ("T790M/WZ4002", "3IKA"), ("T790M", "5GTZ"), ("T790M", "5GNK"),
    ("C797S?", "6JXT"), ("C797S?", "6LUD"), ("C797S?", "5Y9T"),
    ("C797S?", "6JWL"), ("C797S?", "5U8L"), ("L858R", "2ITV"),
    ("L858R/T790M", "6JWL"), ("osi/AZD9291", "4ZAU"),
]

three2one_res = {}
os.makedirs("_probe", exist_ok=True)
parser = PDBParser(QUIET=True)
seen = set()
print(f"{'PDB':6} {'res(A)':7} {'790':4} {'797':4} {'858':4}  ligands")
for state, pid in CANDIDATES:
    if pid in seen:
        continue
    seen.add(pid)
    path = f"_probe/{pid}.pdb"
    try:
        if not os.path.exists(path):
            urllib.request.urlretrieve(f"https://files.rcsb.org/download/{pid}.pdb", path)
    except Exception as e:
        print(f"{pid:6} DOWNLOAD FAIL ({e})")
        continue
    try:
        s = parser.get_structure(pid, path)
        reso = s.header.get("resolution")
        model = s[0]
        # 첫 폴리펩타이드 체인
        chain = None
        for ch in model:
            if any(r.id[0] == " " for r in ch):
                chain = ch; break
        def resname_at(n):
            try:
                r = chain[(" ", n, " ")]
                return r.resname
            except Exception:
                return "-"
        r790, r797, r858 = resname_at(790), resname_at(797), resname_at(858)
        ligs = sorted({r.resname for r in chain if r.id[0].startswith("H_") and r.resname not in ("HOH","NAG","SO4","GOL","EDO","CL","NA","MG","ACT","PO4")})
        print(f"{pid:6} {str(reso):7} {r790:4} {r797:4} {r858:4}  {','.join(ligs)}")
    except Exception as e:
        print(f"{pid:6} PARSE FAIL ({e})")
