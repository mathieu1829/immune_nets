let SessionLoad = 1
let s:so_save = &g:so | let s:siso_save = &g:siso | setg so=0 siso=0 | setl so=-1 siso=-1
let v:this_session=expand("<sfile>:p")
silent only
silent tabonly
cd ~/Documents/python/immune_nets
if expand('%') == '' && !&modified && line('$') <= 1 && getline(1) == ''
  let s:wipebuf = bufnr('%')
endif
let s:shortmess_save = &shortmess
if &shortmess =~ 'A'
  set shortmess=aoOA
else
  set shortmess=aoO
endif
badd +0 src/mappers/networkMapper.py
badd +1 src/analysis/methods/graphStats.py
badd +1 src/models/networkStat.py
badd +1 tests/test_graphStats.py
argglobal
%argdel
edit src/mappers/networkMapper.py
let s:save_splitbelow = &splitbelow
let s:save_splitright = &splitright
set splitbelow splitright
wincmd _ | wincmd |
vsplit
wincmd _ | wincmd |
vsplit
2wincmd h
wincmd w
wincmd w
wincmd _ | wincmd |
split
1wincmd k
wincmd w
let &splitbelow = s:save_splitbelow
let &splitright = s:save_splitright
wincmd t
let s:save_winminheight = &winminheight
let s:save_winminwidth = &winminwidth
set winminheight=0
set winheight=1
set winminwidth=0
set winwidth=1
exe 'vert 1resize ' . ((&columns * 126 + 190) / 381)
exe 'vert 2resize ' . ((&columns * 126 + 190) / 381)
exe '3resize ' . ((&lines * 37 + 38) / 77)
exe 'vert 3resize ' . ((&columns * 127 + 190) / 381)
exe '4resize ' . ((&lines * 37 + 38) / 77)
exe 'vert 4resize ' . ((&columns * 127 + 190) / 381)
argglobal
balt src/mappers/networkMapper.py
setlocal foldmethod=manual
setlocal foldexpr=0
setlocal foldmarker={{{,}}}
setlocal foldignore=#
setlocal foldlevel=0
setlocal foldminlines=1
setlocal foldnestmax=20
setlocal foldenable
silent! normal! zE
let &fdl = &fdl
let s:l = 1 - ((0 * winheight(0) + 37) / 74)
if s:l < 1 | let s:l = 1 | endif
keepjumps exe s:l
normal! zt
keepjumps 1
normal! 0
wincmd w
argglobal
if bufexists(fnamemodify("src/analysis/methods/graphStats.py", ":p")) | buffer src/analysis/methods/graphStats.py | else | edit src/analysis/methods/graphStats.py | endif
if &buftype ==# 'terminal'
  silent file src/analysis/methods/graphStats.py
endif
balt src/mappers/networkMapper.py
setlocal foldmethod=manual
setlocal foldexpr=0
setlocal foldmarker={{{,}}}
setlocal foldignore=#
setlocal foldlevel=0
setlocal foldminlines=1
setlocal foldnestmax=20
setlocal foldenable
silent! normal! zE
let &fdl = &fdl
let s:l = 1 - ((0 * winheight(0) + 37) / 74)
if s:l < 1 | let s:l = 1 | endif
keepjumps exe s:l
normal! zt
keepjumps 1
normal! 0
wincmd w
argglobal
if bufexists(fnamemodify("tests/test_graphStats.py", ":p")) | buffer tests/test_graphStats.py | else | edit tests/test_graphStats.py | endif
if &buftype ==# 'terminal'
  silent file tests/test_graphStats.py
endif
balt src/models/networkStat.py
setlocal foldmethod=manual
setlocal foldexpr=0
setlocal foldmarker={{{,}}}
setlocal foldignore=#
setlocal foldlevel=0
setlocal foldminlines=1
setlocal foldnestmax=20
setlocal foldenable
silent! normal! zE
let &fdl = &fdl
let s:l = 1 - ((0 * winheight(0) + 18) / 36)
if s:l < 1 | let s:l = 1 | endif
keepjumps exe s:l
normal! zt
keepjumps 1
normal! 0
wincmd w
argglobal
if bufexists(fnamemodify("src/models/networkStat.py", ":p")) | buffer src/models/networkStat.py | else | edit src/models/networkStat.py | endif
if &buftype ==# 'terminal'
  silent file src/models/networkStat.py
endif
balt src/analysis/methods/graphStats.py
setlocal foldmethod=manual
setlocal foldexpr=0
setlocal foldmarker={{{,}}}
setlocal foldignore=#
setlocal foldlevel=0
setlocal foldminlines=1
setlocal foldnestmax=20
setlocal foldenable
silent! normal! zE
let &fdl = &fdl
let s:l = 1 - ((0 * winheight(0) + 18) / 36)
if s:l < 1 | let s:l = 1 | endif
keepjumps exe s:l
normal! zt
keepjumps 1
normal! 0
wincmd w
3wincmd w
exe 'vert 1resize ' . ((&columns * 126 + 190) / 381)
exe 'vert 2resize ' . ((&columns * 126 + 190) / 381)
exe '3resize ' . ((&lines * 37 + 38) / 77)
exe 'vert 3resize ' . ((&columns * 127 + 190) / 381)
exe '4resize ' . ((&lines * 37 + 38) / 77)
exe 'vert 4resize ' . ((&columns * 127 + 190) / 381)
tabnext 1
if exists('s:wipebuf') && len(win_findbuf(s:wipebuf)) == 0 && getbufvar(s:wipebuf, '&buftype') isnot# 'terminal'
  silent exe 'bwipe ' . s:wipebuf
endif
unlet! s:wipebuf
set winheight=1 winwidth=20
let &shortmess = s:shortmess_save
let &winminheight = s:save_winminheight
let &winminwidth = s:save_winminwidth
let s:sx = expand("<sfile>:p:r")."x.vim"
if filereadable(s:sx)
  exe "source " . fnameescape(s:sx)
endif
let &g:so = s:so_save | let &g:siso = s:siso_save
set hlsearch
nohlsearch
doautoall SessionLoadPost
unlet SessionLoad
" vim: set ft=vim :
