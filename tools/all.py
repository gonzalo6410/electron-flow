"""Screenshot every compound of the test list (core/compounds.js, typed into the search box as SMILES) and every
reaction scene into contact sheets.
Usage: python3 -I all.py OUTDIR
"""
import asyncio, sys, os, json
from playwright.async_api import async_playwright
from PIL import Image

async def main():
    out = sys.argv[1]; os.makedirs(out, exist_ok=True)
    async with async_playwright() as p:
        b = await p.chromium.launch(executable_path='/usr/bin/google-chrome', args=['--use-gl=swiftshader', '--enable-unsafe-swiftshader', '--no-sandbox'])
        pg = await b.new_page(viewport={'width': 1200, 'height': 800}, device_scale_factor=1)
        errs = []
        pg.on('console', lambda m: errs.append(f'{m.type}: {m.text}') if m.type == 'error' else None)
        pg.on('pageerror', lambda e: errs.append(f'pageerror: {e}'))
        await pg.goto('http://localhost:8765/index.html'); await pg.wait_for_timeout(2500)
        compounds = await pg.evaluate("import('./core/compounds.js').then(m => m.COMPOUNDS.map(c => [c.name, c.smiles]))")
        scenes = await pg.evaluate("[...document.querySelectorAll('#case option')].filter(o => o.value !== '').map(o => [o.value, o.textContent])")
        todo = [('c', name, smi) for name, smi in compounds] + [('s', name, value) for value, name in scenes]
        shots = []
        for q, (kind, name, what) in enumerate(todo):
            if kind == 'c':
                await pg.fill('#search', what); await pg.press('#search', 'Enter')
            else:
                await pg.select_option('#case', what)
            await pg.wait_for_timeout(1300)
            path = f'{out}/{q:03d}.png'
            await pg.locator('#view').screenshot(path=path)
            status = (await pg.inner_text('#status')).replace('\n', ' | ')
            shots.append((name, path, status))
            print(f'{name}: {status}')
        print('\n'.join(e for e in errs if 'favicon' not in e and '404' not in e) or 'no errors')
        await b.close()
    # contact sheets of 12
    for k in range(0, len(shots), 12):
        batch = shots[k:k + 12]
        sheet = Image.new('RGB', (4 * 455, 3 * 400), 'white')
        for q, (name, path, status) in enumerate(batch):
            im = Image.open(path).resize((455, 400))
            sheet.paste(im, ((q % 4) * 455, (q // 4) * 400))
        sheet.save(f'{out}/sheet{k // 12}.png')

asyncio.run(main())
