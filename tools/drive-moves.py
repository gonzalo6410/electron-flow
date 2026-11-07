"""Drive the moves on the Electron Flow page: python3 -I drive-moves.py SCRIPT OUT_PREFIX
Steps (semicolon separated):
  pick:TEXT          choose the picker entry whose name starts with TEXT
  smiles:TEXT        type a name or a SMILES and press enter
  show               click "show the arrows"
  go                 press Enter (make the move)
  button:ID          click a button by id
  wait:MS
  shot:NAME          screenshot and print status, note, error, arrows
  dragpair:I,J       drag from a dot pair of atom I to atom J
  dragbond:I,J,K     drag from the middle of bond I-J to atom K
  dragbondbond:I,J,K,L  drag from the middle of bond I-J to the middle of bond K-L
  key:K
  res:on|off
  eval:JS
"""
import asyncio, sys
from playwright.async_api import async_playwright
PORT = 8765

async def proj(pg, p):
    return await pg.evaluate("""(p) => { const s = window.electronFlow; const v = new s.THREE.Vector3(...p).project(s.camera);
      const r = s.canvas.getBoundingClientRect(); return [r.left + (v.x + 1) / 2 * r.width, r.top + (1 - v.y) / 2 * r.height]; }""", p)
async def atom_xy(pg, i):
    p = await pg.evaluate("(i) => window.electronFlow.positions()[i]", i)
    return await proj(pg, p)
async def pair_xy(pg, i):
    ps = await pg.evaluate("(i) => window.electronFlow.pairs().filter(t => t.index === i).map(t => t.pos)", i)
    if not ps: return None
    return await proj(pg, ps[0])
async def mid_xy(pg, i, j):
    a = await pg.evaluate("(i) => window.electronFlow.positions()[i]", i)
    b = await pg.evaluate("(i) => window.electronFlow.positions()[i]", j)
    return await proj(pg, [(a[0]+b[0])/2, (a[1]+b[1])/2, (a[2]+b[2])/2])
async def drag(pg, a, b):
    await pg.mouse.move(a[0], a[1]); await pg.mouse.down()
    steps = 12
    for q in range(1, steps + 1):
        await pg.mouse.move(a[0] + (b[0]-a[0]) * q / steps, a[1] + (b[1]-a[1]) * q / steps); await pg.wait_for_timeout(30)
    await pg.mouse.up(); await pg.wait_for_timeout(300)

async def main():
    script, prefix = sys.argv[1], sys.argv[2]
    async with async_playwright() as p:
        b = await p.chromium.launch(executable_path='/usr/bin/google-chrome', args=['--use-gl=swiftshader', '--enable-unsafe-swiftshader', '--no-sandbox'])
        pg = await b.new_page(viewport={'width': 1200, 'height': 800}, device_scale_factor=2)
        msgs = []
        pg.on('console', lambda m: msgs.append(f'{m.type}: {m.text}') if m.type == 'error' else None)
        pg.on('pageerror', lambda e: msgs.append(f'pageerror: {e}'))
        await pg.goto(f'http://localhost:{PORT}/index.html'); await pg.wait_for_timeout(2500)
        for step in script.split(';'):
            step = step.strip()
            if not step: continue
            cmd, _, arg = step.partition(':')
            if cmd == 'pick':
                await pg.evaluate("""(t) => { const sel = document.getElementById('case'); const o = [...sel.options].find(o => o.textContent.startsWith(t)); sel.value = o.value; sel.dispatchEvent(new Event('change')); }""", arg)
                await pg.wait_for_timeout(1500)
            elif cmd == 'smiles':
                await pg.fill('#search', arg); await pg.press('#search', 'Enter'); await pg.wait_for_timeout(1200)
            elif cmd == 'show': await pg.click('#show'); await pg.wait_for_timeout(400)
            elif cmd == 'addc':
                await pg.evaluate("""(t) => { const b = document.getElementById('search'); b.focus(); b.value = t; b.dispatchEvent(new Event('input')); b.dispatchEvent(new KeyboardEvent('keydown', { key: 'Enter', shiftKey: true, bubbles: true })); }""", arg)
                await pg.wait_for_timeout(1500)
            elif cmd == 'go': await pg.keyboard.press('Enter'); await pg.wait_for_timeout(100)
            elif cmd == 'button': await pg.click('#' + arg); await pg.wait_for_timeout(400)
            elif cmd == 'wait': await pg.wait_for_timeout(int(arg))
            elif cmd == 'key': await pg.keyboard.press(arg); await pg.wait_for_timeout(400)
            elif cmd == 'res':
                checked = await pg.is_checked('#resonance')
                if (arg == 'on') != checked: await pg.click('.switch')
                await pg.wait_for_timeout(300)
            elif cmd == 'dragpair':
                i, j = [int(x) for x in arg.split(',')]
                a = await pair_xy(pg, i); bb = await atom_xy(pg, j)
                if not a: print('no pair on', i); continue
                await drag(pg, a, bb)
            elif cmd == 'dragbond':
                i, j, k = [int(x) for x in arg.split(',')]
                await drag(pg, await mid_xy(pg, i, j), await atom_xy(pg, k))
            elif cmd == 'dragbondbond':
                i, j, k, l = [int(x) for x in arg.split(',')]
                await drag(pg, await mid_xy(pg, i, j), await mid_xy(pg, k, l))
            elif cmd == 'eval':
                print('eval:', await pg.evaluate(arg))
            elif cmd == 'shot':
                await pg.screenshot(path=f'{prefix}-{arg}.png')
                status = (await pg.inner_text('#status')).replace('\n', ' | ')
                note = await pg.inner_text('#note'); err = await pg.inner_text('#error')
                arrows = await pg.evaluate("JSON.stringify(window.electronFlow.arrows())")
                print(f'[{arg}] status: {status}\n[{arg}] note: {note}\n[{arg}] arrows: {arrows}' + (f'\n[{arg}] error: {err}' if err else ''))
        if msgs: print('CONSOLE:', *msgs, sep='\n  ')
        await b.close()
asyncio.run(main())
