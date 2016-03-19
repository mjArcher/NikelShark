let SessionLoad = 1
if &cp | set nocp | endif
let s:so_save = &so | let s:siso_save = &siso | set so=0 siso=0
let v:this_session=expand("<sfile>:p")
silent only
cd ~/Code/computing/solid/solid-1D/src
if expand('%') == '' && !&modified && line('$') <= 1 && getline(1) == ''
  let s:wipebuf = bufnr('%')
endif
set shortmess=aoO
badd +591 Elastic1D.cpp
badd +0 Input.hpp
badd +0 InputSolid.h
badd +45 ~/Code/computing/solid/solid-1D/settings/barton1D.cfg
badd +0 ~/Code/computing/solid/solid-1D/settings/input.cfg
badd +0 Input.cpp
badd +0 Makefile
badd +51 ~/Code/computing/solid/solid-1D/settings/simple.cfg
args Elastic1D.cpp
edit Elastic1D.cpp
set splitbelow splitright
wincmd _ | wincmd |
vsplit
wincmd _ | wincmd |
vsplit
2wincmd h
wincmd w
wincmd w
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
exe 'vert 1resize ' . ((&columns * 79 + 119) / 239)
exe 'vert 2resize ' . ((&columns * 79 + 119) / 239)
exe 'vert 3resize ' . ((&columns * 79 + 119) / 239)
argglobal
let s:l = 572 - ((38 * winheight(0) + 31) / 62)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
572
normal! 02l
wincmd w
argglobal
edit InputSolid.h
let s:l = 151 - ((30 * winheight(0) + 31) / 62)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
151
normal! 034l
wincmd w
argglobal
edit ~/Code/computing/solid/solid-1D/settings/barton1D.cfg
let s:l = 22 - ((21 * winheight(0) + 31) / 62)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
22
normal! 0
wincmd w
2wincmd w
exe 'vert 1resize ' . ((&columns * 79 + 119) / 239)
exe 'vert 2resize ' . ((&columns * 79 + 119) / 239)
exe 'vert 3resize ' . ((&columns * 79 + 119) / 239)
tabedit Input.cpp
set splitbelow splitright
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
argglobal
let s:l = 108 - ((43 * winheight(0) + 31) / 62)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
108
normal! 01l
2wincmd w
tabedit ~/Code/computing/solid/solid-1D/settings/input.cfg
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
exe 'vert 1resize ' . ((&columns * 119 + 119) / 239)
exe 'vert 2resize ' . ((&columns * 119 + 119) / 239)
argglobal
let s:l = 12 - ((11 * winheight(0) + 31) / 62)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
12
normal! 0
wincmd w
argglobal
edit Input.hpp
let s:l = 259 - ((33 * winheight(0) + 31) / 62)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
259
normal! 028l
wincmd w
2wincmd w
exe 'vert 1resize ' . ((&columns * 119 + 119) / 239)
exe 'vert 2resize ' . ((&columns * 119 + 119) / 239)
tabedit Makefile
set splitbelow splitright
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
argglobal
let s:l = 25 - ((24 * winheight(0) + 31) / 62)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
25
normal! 038l
2wincmd w
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
