#!/usr/bin/env node
// Fetch APS papers via headed Playwright browser (TIB VPN required).
// Usage: node scripts/fetch_aps_papers.mjs
//
// 1. Ensure TIB VPN is active.
// 2. Script opens headed Chromium and navigates to an APS page.
// 3. Click the Cloudflare challenge in the browser window.
// 4. Script auto-detects when you pass it and fetches all 9 PDFs.

import { chromium } from 'playwright';
import { writeFileSync, existsSync, mkdirSync } from 'fs';
import { resolve } from 'path';

const PAPERS_DIR = resolve(import.meta.dirname, '..', 'sources', 'papers');

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

async function waitForCloudflare(page, timeoutMs = 180000) {
  console.log('Waiting for Cloudflare challenge to be solved...');
  console.log('>>> Click the challenge in the browser window <<<');
  console.log(`>>> You have ${timeoutMs/1000}s — take your time <<<\n`);

  // Wait for an element that only exists on the real APS page, not on Cloudflare.
  // APS abstract pages have <h3 class="title"> or <div id="article"> etc.
  // We try multiple selectors and also check for any <meta> tag with "aps.org".
  try {
    await page.waitForFunction(() => {
      // Check for real APS content — any of these signals we're past Cloudflare
      return document.querySelector('meta[name="citation_title"]')
        || document.querySelector('h3.title')
        || document.querySelector('#article')
        || document.querySelector('.article-title')
        || (document.title && document.title.includes('Phys. Rev.'));
    }, { timeout: timeoutMs });

    const title = await page.title();
    console.log(`Page loaded: "${title}"`);
    console.log('Cloudflare passed!\n');
    await new Promise(r => setTimeout(r, 2000));
    return true;
  } catch (_) {
    return false;
  }
}

async function main() {
  mkdirSync(PAPERS_DIR, { recursive: true });

  console.log('Launching HEADED Chromium (persistent profile)...');
  console.log('Make sure TIB VPN is active!\n');

  // Use persistent context — Cloudflare cookies survive across runs
  const userDataDir = resolve(PAPERS_DIR, '..', '.browser-profile');
  mkdirSync(userDataDir, { recursive: true });

  const context = await chromium.launchPersistentContext(userDataDir, {
    headless: false,
    args: ['--disable-blink-features=AutomationControlled'],
    viewport: { width: 1280, height: 900 },
  });
  const page = context.pages()[0] || await context.newPage();

  // Navigate to an APS abstract page to trigger Cloudflare
  const triggerUrl = 'https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.51.605';
  console.log(`Navigating to: ${triggerUrl}`);
  await page.goto(triggerUrl, { waitUntil: 'domcontentloaded', timeout: 60000 });

  const passed = await waitForCloudflare(page);
  if (!passed) {
    console.error('Timed out waiting for Cloudflare. Aborting.');
    await context.close();
    process.exit(1);
  }

  console.log('Fetching PDFs...\n');

  let downloaded = 0;
  let failed = 0;
  let skipped = 0;

  for (const paper of PAPERS) {
    const outPath = resolve(PAPERS_DIR, paper.file);
    if (existsSync(outPath)) {
      console.log(`SKIP ${paper.id}: ${paper.file} (already exists)`);
      skipped++;
      continue;
    }

    try {
      process.stdout.write(`FETCH ${paper.id}: ${paper.file} ... `);
      const response = await page.request.get(paper.url, { timeout: 30000 });

      if (response.status() !== 200) {
        console.log(`FAIL (HTTP ${response.status()})`);
        failed++;
        continue;
      }

      const body = await response.body();
      const header = body.slice(0, 5).toString();
      if (header !== '%PDF-') {
        console.log(`FAIL (not a PDF, got: ${header})`);
        failed++;
        continue;
      }

      writeFileSync(outPath, body);
      console.log(`OK (${(body.length/1024).toFixed(0)} KB)`);
      downloaded++;

      await new Promise(r => setTimeout(r, 1500));
    } catch (e) {
      console.log(`ERROR: ${e.message}`);
      failed++;
    }
  }

  console.log(`\nDone: ${downloaded} downloaded, ${failed} failed, ${skipped} skipped`);
  await context.close();
}

main().catch(e => { console.error(e); process.exit(1); });
