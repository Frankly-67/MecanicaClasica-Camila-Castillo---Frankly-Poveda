import sys, json, base64, subprocess, os, asyncio
from playwright.async_api import async_playwright
from cfg import build
FPS=30
async def main(name, out, only=None):
    cfg=build(); beats=cfg[name]['beats']
    frames=[]  # (beat,tt,tAbs)
    t=0
    for bi,b in enumerate(beats):
        n=int(round(b['dur']*FPS))
        for k in range(n): frames.append((bi,k/FPS,t+k/FPS))
        t+=b['dur']
    if only is not None:
        sel=[]
        for spec in only:   # "beat:frac"
            bi,fr=spec.split(':'); bi=int(bi); fr=float(fr)
            b=beats[bi]; tt=b['dur']*fr; sel.append((bi,tt,tt))
        frames=sel
    async with async_playwright() as p:
        br=await p.chromium.launch(); pg=await br.new_page(viewport={'width':1280,'height':720})
        await pg.goto('file:///home/claude/video-anim/escena.html')
        await pg.evaluate("c=>{window.CFG=c}",cfg)
        errs=[]; pg.on('pageerror',lambda e:errs.append(str(e))); pg.on('console',lambda m: errs.append(m.text) if m.type=='error' else None)
        if only is not None:
            for i,(bi,tt,ta) in enumerate(frames):
                await pg.evaluate("a=>renderFrame(a[0],CFG,a[1],a[2],a[3])",[name,bi,tt,ta])
                d=await pg.evaluate("document.getElementById('c').toDataURL('image/jpeg',0.92)")
                open(f'{out}_{i}.jpg','wb').write(base64.b64decode(d.split(',')[1]))
            print('errs',errs[:5]); await br.close(); return
        ff=subprocess.Popen(['ffmpeg','-y','-loglevel','error','-f','image2pipe','-framerate',str(FPS),'-c:v','mjpeg','-i','-','-c:v','libx264','-crf','19','-preset','veryfast','-pix_fmt','yuv420p',out+'.tmp.mp4'],stdin=subprocess.PIPE)
        for i,(bi,tt,ta) in enumerate(frames):
            await pg.evaluate("a=>renderFrame(a[0],CFG,a[1],a[2],a[3])",[name,bi,tt,ta])
            d=await pg.evaluate("document.getElementById('c').toDataURL('image/jpeg',0.92)")
            ff.stdin.write(base64.b64decode(d.split(',')[1]))
        ff.stdin.close(); rc=ff.wait()
        if rc==0: os.rename(out+'.tmp.mp4',out)
        print(name,'rc',rc,'errs',errs[:5]); await br.close()
if __name__=='__main__':
    name=sys.argv[1]; out=sys.argv[2]; only=sys.argv[3:] or None
    asyncio.run(main(name,out,only))
