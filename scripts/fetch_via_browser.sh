#!/bin/bash
# Download APS papers via the existing playwright-cli headed browser session.
# Requires: browser already open and authenticated (Cloudflare passed, TIB VPN active).

PAPERS_DIR="/home/tobiasosborne/Projects/FQHE/sources/papers"

declare -A PAPERS
PAPERS[P01_Haldane_PRL_51_605_1983.pdf]="https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.51.605"
PAPERS[P02_Laughlin_PRL_50_1395_1983.pdf]="https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.50.1395"
PAPERS[P03_Jain_PRL_63_199_1989.pdf]="https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.63.199"
PAPERS[P04_Fano_PRB_34_2670_1986.pdf]="https://journals.aps.org/prb/pdf/10.1103/PhysRevB.34.2670"
PAPERS[P06_Tsui_PRL_48_1559_1982.pdf]="https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.48.1559"
PAPERS[P07_Willett_PRL_59_1776_1987.pdf]="https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.59.1776"
PAPERS[P10_Dykhne_PRB_50_2369_1994.pdf]="https://journals.aps.org/prb/pdf/10.1103/PhysRevB.50.2369"
PAPERS[P12_Stormer_RMP_71_875_1999.pdf]="https://journals.aps.org/rmp/pdf/10.1103/RevModPhys.71.875"
PAPERS[P15_Haldane_Rezayi_PRL_54_237_1985.pdf]="https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.54.237"

OK=0
FAIL=0
SKIP=0

for FNAME in $(echo "${!PAPERS[@]}" | tr ' ' '\n' | sort); do
  URL="${PAPERS[$FNAME]}"
  OUTPATH="${PAPERS_DIR}/${FNAME}"

  if [ -f "$OUTPATH" ]; then
    SIZE=$(stat -c%s "$OUTPATH")
    FTYPE=$(file -b "$OUTPATH" | head -c 3)
    if [ "$SIZE" -gt 10000 ] && [ "$FTYPE" = "PDF" ]; then
      echo "SKIP $FNAME (${SIZE} bytes, valid PDF)"
      SKIP=$((SKIP+1))
      continue
    else
      echo "REPLACING $FNAME (invalid: ${SIZE} bytes, type: $FTYPE)"
      rm "$OUTPATH"
    fi
  fi

  echo -n "FETCH $FNAME ... "

  RAW=$(playwright-cli run-code "async page => {
    const resp = await page.request.get('${URL}', { timeout: 30000 });
    if (resp.status() !== 200) return 'ERROR:' + resp.status();
    const body = await resp.body();
    return body.toString('base64');
  }" 2>&1)

  B64=$(echo "$RAW" | grep -A1 '### Result' | tail -1 | tr -d '"')

  if [[ "$B64" == ERROR* ]] || [[ -z "$B64" ]] || [[ "$B64" == *"### Error"* ]]; then
    echo "FAIL ($B64)"
    FAIL=$((FAIL+1))
    continue
  fi

  echo "$B64" | base64 -d > "$OUTPATH"
  SIZE=$(stat -c%s "$OUTPATH")
  FTYPE=$(file -b "$OUTPATH" | head -c 20)

  if [[ "$FTYPE" == PDF* ]] && [ "$SIZE" -gt 10000 ]; then
    echo "OK (${SIZE} bytes, ${FTYPE})"
    OK=$((OK+1))
  else
    echo "FAIL (bad file: ${SIZE} bytes, ${FTYPE})"
    rm "$OUTPATH"
    FAIL=$((FAIL+1))
  fi

  sleep 2
done

echo ""
echo "Done: ${OK} downloaded, ${FAIL} failed, ${SKIP} skipped"
