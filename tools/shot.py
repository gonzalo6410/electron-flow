# Electron Flow. Copyright (c) 2026 Gonzalo Fernandez de la Mora. All rights reserved (see LICENSE).
import asyncio, sys
from playwright.async_api import async_playwright
async def main():
    async with async_playwright() as p:
        b = await p.chromium.launch(executable_path='/usr/bin/google-chrome', args=['--use-gl=swiftshader','--enable-unsafe-swiftshader','--no-sandbox'])
        pg = await b.new_page(viewport={'width':1200,'height':800}, device_scale_factor=2)
        msgs=[]
        pg.on('console', lambda m: msgs.append(f'{m.type}: {m.text}') if m.type=='error' else None)
        pg.on('pageerror', lambda e: msgs.append(f'pageerror: {e}'))
        await pg.goto('http://localhost:8765/index.html')
        await pg.wait_for_timeout(2500)
        smi = sys.argv[1]; res = sys.argv[2] == 'on'; wait = int(sys.argv[3]); out = sys.argv[4]
        await pg.fill('#search', smi); await pg.press('#search', 'Enter')
        if res: await pg.click('.switch')
        await pg.wait_for_timeout(wait)
        await pg.screenshot(path=out)
        print('\n'.join(m for m in msgs if 'favicon' not in m and '404' not in m) or 'no errors')
        print(await pg.inner_text('#status'))
        await b.close()
asyncio.run(main())
