#!/usr/bin/env node
// Fetch APS papers using Playwright with saved browser cookies (TIB VPN auth)
// Usage: node scripts/fetch_aps_papers.mjs

import { chromium } from '/home/tobiasosborne/.npm/_npx/e41f203b7505f1fb/node_modules/playwright-core/index.mjs';
import { readFileSync, writeFileSync, existsSync } from 'fs';
import { resolve } from 'path';

const PAPERS_DIR = resolve(import.meta.dirname, '..', 'sources', 'papers');
const STATE_FILE = resolve(import.meta.dirname, '..', 'sources', '.aps-cookies.json');

const PAPERS = [
  { id: 'P01', file: 'P01_Haldane_PRL_51_605_1983.pdf',       url: 'https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.51.605' },
  { id: 'P02', file: 'P02_Laughlin_PRL_50_1395_1983.pdf',      url: 'https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.50.1395' },
  { id: 'P03', file: 'P03_Jain_PRL_63_199_1989.pdf',           url: 'https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.63.199' },
  { id: 'P04', file: 'P04_Fano_PRB_34_2670_1986.pdf',          url: 'https://journals.aps.org/prb/pdf/10.1103/PhysRevB.34.2670' },
  { id: 'P06', file: 'P06_Tsui_PRL_48_1559_1982.pdf',          url: 'https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.48.1559' },
  { id: 'P07', file: 'P07_Willett_PRL_59_1776_1987.pdf',       url: 'https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.59.1776' },
  { id: 'P10', file: 'P10_Dykhne_PRB_50_2369_1994.pdf',        url: 'https://journals.aps.org/prb/pdf/10.1103/PhysRevB.50.2369' },
  { id: 'P12', file: 'P12_Stormer_RMP_71_875_1999.pdf',        url: 'https://journals.aps.org/rmp/pdf/10.1103/RevModPhys.71.875' },
  { id: 'P15', file: 'P15_Haldane_Rezayi_PRL_54_237_1985.pdf', url: 'https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.54.237' },
];

async function main() {
  if (!existsSync(STATE_FILE)) {
    console.error('No cookie state file found at', STATE_FILE);
    console.error('Run playwright-cli in headed mode first to pass Cloudflare, then state-save.');
    process.exit(1);
  }

  // Find chromium
  const browserPaths = [
    '/home/tobiasosborne/.cache/ms-playwright/chromium-1212/chrome-linux64/chrome',
    '/home/tobiasosborne/.cache/ms-playwright/chromium-1208/chrome-linux64/chrome',
  ];
  let execPath = browserPaths.find(p => existsSync(p));
  if (!execPath) {
    console.error('No chromium found');
    process.exit(1);
  }

  console.log('Launching browser...');
  const browser = await chromium.launch({
    executablePath: execPath,
    headless: true
  });

  const context = await browser.newContext({
    storageState: STATE_FILE,
    userAgent: 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/136.0.0.0 Safari/537.36',
  });

  const page = await context.newPage();

  let downloaded = 0;
  let failed = 0;

  for (const paper of PAPERS) {
    const outPath = resolve(PAPERS_DIR, paper.file);
    if (existsSync(outPath)) {
      console.log(`SKIP ${paper.id}: ${paper.file} (already exists)`);
      continue;
    }

    try {
      console.log(`FETCH ${paper.id}: ${paper.url}`);
      const response = await page.request.get(paper.url, { timeout: 30000 });

      if (response.status() !== 200) {
        console.error(`  FAIL ${paper.id}: HTTP ${response.status()}`);
        failed++;
        continue;
      }

      const body = await response.body();

      // Verify it's actually a PDF
      const header = body.slice(0, 5).toString();
      if (header !== '%PDF-') {
        console.error(`  FAIL ${paper.id}: Not a PDF (got: ${header})`);
        failed++;
        continue;
      }

      writeFileSync(outPath, body);
      console.log(`  OK ${paper.id}: ${paper.file} (${(body.length/1024).toFixed(0)} KB)`);
      downloaded++;

      // Small delay to be polite
      await new Promise(r => setTimeout(r, 1500));

    } catch (e) {
      console.error(`  FAIL ${paper.id}: ${e.message}`);
      failed++;
    }
  }

  await browser.close();
  console.log(`\nDone: ${downloaded} downloaded, ${failed} failed, ${PAPERS.length - downloaded - failed} skipped`);
}

main().catch(e => { console.error(e); process.exit(1); });
