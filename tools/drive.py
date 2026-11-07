# Electron Flow. Copyright (c) 2026 Gonzalo Fernandez de la Mora. All rights reserved (see LICENSE).
"""Drive the Electron Flow page with playwright and take screenshots.

Usage: python3 -I drive.py SCRIPT OUT_PREFIX
SCRIPT is a semicolon-separated list of steps:
  smiles:TEXT        type a name or a SMILES and press enter
  res:on|off         set the resonance switch
  wait:MS            wait
  shot:NAME          save OUT_PREFIX-NAME.png and print the status and note
  click:I            click atom I (by its index) in the picture
  shiftclick:I       shift-click atom I
  key:K              press a key on the page (e.g. c, Delete, +, Shift+H, ArrowRight, Space)
  drag:DX,DY         drag the picture by DX,DY pixels from its centre
  fill               click the fill hydrogens button
  clear              click the new molecule button
  table:EL           click element EL in the periodic table
  drop:EL,I          drag element EL from the table onto atom I (I = -1 for empty space)
  hover:I            move the pointer over atom I and print the tip
  eval:JS            evaluate JS in the page and print the result
"""
import asyncio, sys, json
from playwright.async_api import async_playwright

PORT = 8765

async def atom_xy(pg, i):
    # project atom i to screen coordinates through the page's own camera
    return await pg.evaluate("""(i) => {
      const s = window.electronFlow; if (!s) return null;
      const p = s.positions()[i]; if (!p) return null;
      const v = new s.THREE.Vector3(...p).project(s.camera);
      const r = s.canvas.getBoundingClientRect();
      return [r.left + (v.x + 1) / 2 * r.width, r.top + (1 - v.y) / 2 * r.height];
    }""", i)

async def main():
    script, prefix = sys.argv[1], sys.argv[2]
    async with async_playwright() as p:
        b = await p.chromium.launch(executable_path='/usr/bin/google-chrome', args=['--use-gl=swiftshader', '--enable-unsafe-swiftshader', '--no-sandbox'])
        pg = await b.new_page(viewport={'width': 1200, 'height': 800}, device_scale_factor=2)
        msgs = []
        pg.on('console', lambda m: msgs.append(f'{m.type}: {m.text}') if m.type == 'error' else None)
        pg.on('pageerror', lambda e: msgs.append(f'pageerror: {e}'))
        await pg.goto(f'http://localhost:{PORT}/index.html')
        await pg.wait_for_timeout(2500)
        for step in script.split(';'):
            step = step.strip()
            if not step: continue
            cmd, _, arg = step.partition(':')
            if cmd == 'smiles':
                await pg.fill('#search', arg); await pg.press('#search', 'Enter'); await pg.wait_for_timeout(1200)
            elif cmd == 'res':
                checked = await pg.is_checked('#resonance')
                if (arg == 'on') != checked: await pg.click('.switch')
                await pg.wait_for_timeout(300)
            elif cmd == 'wait':
                await pg.wait_for_timeout(int(arg))
            elif cmd == 'shot':
                await pg.screenshot(path=f'{prefix}-{arg}.png')
                status = (await pg.inner_text('#status')).replace('\n', ' | ')
                note = await pg.inner_text('#note')
                err = await pg.inner_text('#error')
                print(f'[{arg}] status: {status}')
                print(f'[{arg}] note: {note}')
                if err: print(f'[{arg}] error: {err}')
            elif cmd in ('click', 'shiftclick'):
                xy = await atom_xy(pg, int(arg))
                if not xy: print(f'no atom {arg}'); continue
                mods = ['Shift'] if cmd == 'shiftclick' else []
                await pg.mouse.click(xy[0], xy[1], modifiers=mods) if False else None
                if mods: await pg.keyboard.down('Shift')
                await pg.mouse.move(xy[0], xy[1]); await pg.mouse.down(); await pg.mouse.up()
                if mods: await pg.keyboard.up('Shift')
                await pg.wait_for_timeout(300)
            elif cmd == 'key':
                await pg.keyboard.press(arg); await pg.wait_for_timeout(500)
            elif cmd == 'drag':
                dx, dy = [int(x) for x in arg.split(',')]
                box = await pg.locator('#view canvas').bounding_box()
                cx, cy = box['x'] + box['width'] / 2, box['y'] + box['height'] / 2
                await pg.mouse.move(cx, cy); await pg.mouse.down()
                for k in range(1, 11):
                    await pg.mouse.move(cx + dx * k / 10, cy + dy * k / 10)
                await pg.mouse.up(); await pg.wait_for_timeout(800)
            elif cmd == 'fill':
                await pg.click('#fill'); await pg.wait_for_timeout(500)
            elif cmd == 'clear':
                await pg.click('#clear'); await pg.wait_for_timeout(500)
            elif cmd == 'table':
                await pg.click(f'#table button:text-is("{arg}")'); await pg.wait_for_timeout(500)
            elif cmd == 'drop':
                el, i = arg.split(','); i = int(i)
                src = await pg.locator(f'#table button:text-is("{el}")').bounding_box()
                if i >= 0: xy = await atom_xy(pg, i)
                else:
                    box = await pg.locator('#view canvas').bounding_box(); xy = [box['x'] + 60, box['y'] + 60]
                # synthesise the drag events the page listens for
                await pg.evaluate("""([el, x, y]) => {
                  const dt = new DataTransfer(); dt.setData('text/element', el);
                  const c = document.querySelector('#view canvas');
                  c.dispatchEvent(new DragEvent('dragover', { bubbles: true, clientX: x, clientY: y, dataTransfer: dt }));
                  c.dispatchEvent(new DragEvent('drop', { bubbles: true, clientX: x, clientY: y, dataTransfer: dt }));
                }""", [el, xy[0], xy[1]])
                await pg.wait_for_timeout(500)
            elif cmd == 'hover':
                xy = await atom_xy(pg, int(arg))
                if not xy: print(f'no atom {arg}'); continue
                await pg.mouse.move(xy[0], xy[1]); await pg.wait_for_timeout(300)
                print(f'[hover {arg}] ' + (await pg.inner_text('#tip')).replace('\n', ' / '))
            elif cmd == 'eval':
                print('[eval] ' + json.dumps(await pg.evaluate(arg)))
        errs = [m for m in msgs if 'favicon' not in m and '404' not in m]
        print('\n'.join(errs) or 'no errors')
        await b.close()

asyncio.run(main())
