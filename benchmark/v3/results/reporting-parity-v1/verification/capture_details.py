from pathlib import Path
import json
from playwright.sync_api import sync_playwright
root=Path(__file__).parent/'replay-v4'
with sync_playwright() as p:
 browser=p.chromium.launch(headless=True);page=browser.new_page(viewport={'width':1440,'height':1000})
 folder=root/f'{5:032x}'
 page.goto((folder/'gspa_out/sample_browser.html').as_uri());page.locator('[data-tab="ledger"]').click()
 row=page.locator('#ledger-body table').first.locator('tr').filter(has_text='Metabolic experiments')
 row.locator('summary').click()
 row.locator('pre').filter(has_text='"status": "UNKNOWN"').first.scroll_into_view_if_needed()
 page.screenshot(path=str(folder/'unknown-details.png'))
 page.close()
 for i in [1,2]:
  folder=root/f'{i:032x}';page=browser.new_page(viewport={'width':1440,'height':1000})
  data=json.loads((folder/'status.json').read_text());data['ledger']=json.loads((folder/'summary-after.json').read_text())
  page.route('http://gspa.test/**',lambda route:route.fulfill(status=200,body=(Path('/home/leechuck/Public/software/gspa/service/web/index.html').read_text() if route.request.url=='http://gspa.test/' else json.dumps(data)),content_type='text/html' if route.request.url=='http://gspa.test/' else 'application/json'))
  page.goto('http://gspa.test/#'+f'{i:032x}')
  page.wait_for_function("document.querySelector('#ledger').textContent.includes('Scoring:')")
  page.locator('#ledger').scroll_into_view_if_needed();page.screenshot(path=str(folder/'portal-details.png'));page.close()
 browser.close()
