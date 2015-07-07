let SessionLoad = 1
if &cp | set nocp | endif
let s:so_save = &so | let s:siso_save = &siso | set so=0 siso=0
let v:this_session=expand("<sfile>:p")
silent only
cd ~/Dropbox/2013-2014/Code/Solid/myCode/dev
if expand('%') == '' && !&modified && line('$') <= 1 && getline(1) == ''
  let s:wipebuf = bufnr('%')
endif
set shortmess=aoO
badd +32 /local/data/public/ma595/Old/MiniProjects/MiniProject1/CODE/EulerSolver/a_slic_res.sh
badd +65 ~/Dropbox/2013-2014/Code/Solid/myCode/SimpleElastic/Scripts/Solid1DTestCase.sh
badd +1 ~/Dropbox/2013-2014/Code/Solid/myCode/SimpleElastic/Scripts/test_case.sh
badd +1 ~/Dropbox/2013-2014/Code/Solid/myCode/SimpleElastic/Scripts/examples/latexplot.plt
badd +0 ~/Dropbox/2013-2014/Code/Solid/myCode/SimpleElastic/Scripts/examples/comparePlots.sh
badd +0 ~/Dropbox/2013-2014/Code/Solid/myCode/SimpleElastic/Scripts/examples/Plot.sh
badd +0 ~/Dropbox/2013-2014/Code/Solid/myCode/SimpleElastic/Scripts/examples/gnuplotExample.sh
badd +1 ~/Dropbox/2013-2014/Code/Solid/myCode/SimpleElastic/Scripts/latexplot.plt
badd +28 ~/Dropbox/2013-2014/Code/Solid/myCode/SimpleElastic/Scripts/profile.plt
badd +8 ~/Dropbox/2013-2014/Code/Solid/myCode/SimpleElastic/Scripts/profileLatex.plt
badd +0 ~/Dropbox/2013-2014/Code/Solid/myCode/SimpleElastic/Scripts/Solid1DAll.sh
badd +1 ~/Dropbox/2013-2014/Code/Solid/myCode/SimpleElastic/Scripts/profile
badd +9 ./Scripts/Solid1DAll.sh
badd +0 ./Scripts/Solid1DTestCase.sh
args /local/data/public/ma595/Old/MiniProjects/MiniProject1/CODE/EulerSolver/a_slic_res.sh
edit ./Scripts/Solid1DTestCase.sh
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
exe 'vert 1resize ' . ((&columns * 104 + 104) / 209)
exe 'vert 2resize ' . ((&columns * 104 + 104) / 209)
argglobal
let s:l = 40 - ((36 * winheight(0) + 27) / 55)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
40
normal! 019l
wincmd w
argglobal
edit ./Scripts/Solid1DTestCase.sh
let s:l = 189 - ((53 * winheight(0) + 27) / 55)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
189
normal! 0
wincmd w
exe 'vert 1resize ' . ((&columns * 104 + 104) / 209)
exe 'vert 2resize ' . ((&columns * 104 + 104) / 209)
tabnew
set splitbelow splitright
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
argglobal
enew
tabedit ~/Dropbox/2013-2014/Code/Solid/myCode/SimpleElastic/Scripts/Solid1DAll.sh
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
exe 'vert 1resize ' . ((&columns * 100 + 104) / 209)
exe 'vert 2resize ' . ((&columns * 108 + 104) / 209)
argglobal
let s:l = 9 - ((8 * winheight(0) + 27) / 55)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
9
normal! 029l
wincmd w
argglobal
edit ~/Dropbox/2013-2014/Code/Solid/myCode/SimpleElastic/Scripts/profile.plt
let s:l = 1 - ((0 * winheight(0) + 27) / 55)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
1
normal! 0
wincmd w
exe 'vert 1resize ' . ((&columns * 100 + 104) / 209)
exe 'vert 2resize ' . ((&columns * 108 + 104) / 209)
tabedit ~/Dropbox/2013-2014/Code/Solid/myCode/SimpleElastic/Scripts/examples/Plot.sh
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
exe 'vert 1resize ' . ((&columns * 100 + 104) / 209)
exe 'vert 2resize ' . ((&columns * 108 + 104) / 209)
argglobal
let s:l = 8 - ((7 * winheight(0) + 27) / 55)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
8
normal! 0
wincmd w
argglobal
edit ~/Dropbox/2013-2014/Code/Solid/myCode/SimpleElastic/Scripts/examples/comparePlots.sh
let s:l = 67 - ((54 * winheight(0) + 27) / 55)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
67
normal! 0
wincmd w
exe 'vert 1resize ' . ((&columns * 100 + 104) / 209)
exe 'vert 2resize ' . ((&columns * 108 + 104) / 209)
tabedit ~/Dropbox/2013-2014/Code/Solid/myCode/SimpleElastic/Scripts/examples/latexplot.plt
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
exe 'vert 1resize ' . ((&columns * 100 + 104) / 209)
exe 'vert 2resize ' . ((&columns * 108 + 104) / 209)
argglobal
let s:l = 8 - ((7 * winheight(0) + 27) / 55)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
8
normal! 0
wincmd w
argglobal
edit ~/Dropbox/2013-2014/Code/Solid/myCode/SimpleElastic/Scripts/examples/gnuplotExample.sh
let s:l = 54 - ((45 * winheight(0) + 27) / 55)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
54
normal! 0
wincmd w
exe 'vert 1resize ' . ((&columns * 100 + 104) / 209)
exe 'vert 2resize ' . ((&columns * 108 + 104) / 209)
tabnext 1
if exists('s:wipebuf')
  silent exe 'bwipe ' . s:wipebuf
endif
unlet! s:wipebuf
set winheight=1 winwidth=20 shortmess=filnxtToOI
let s:sx = expand("<sfile>:p:r")."x.vim"
if file_readable(s:sx)
  exe "source " . fnameescape(s:sx)
endif
let &so = s:so_save | let &siso = s:siso_save
doautoall SessionLoadPost
unlet SessionLoad
" vim: set ft=vim :
