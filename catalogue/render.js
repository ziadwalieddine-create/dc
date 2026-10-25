const { chromium } = require('playwright');
(async () => {
  const url = process.argv[2];
  const proxy = process.env.HTTPS_PROXY || process.env.https_proxy;
  const browser = await chromium.launch({ proxy: proxy ? { server: proxy } : undefined });
  const page = await browser.newPage();
  try {
    try { await page.goto(url, { waitUntil: 'networkidle', timeout: 30000 }); }
    catch (e) { await page.goto(url, { waitUntil: 'domcontentloaded', timeout: 45000 }); await page.waitForTimeout(6000); }
    await page.waitForTimeout(1500);
    const txt = await page.evaluate(() => document.body.innerText);
    console.log(txt.slice(0, +process.argv[3] || 3000));
  } catch (e) { console.log('ERR', e.message.slice(0, 300)); }
  await browser.close();
})();
