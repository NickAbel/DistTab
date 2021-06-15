R = [4,4];
d = ones(size(R));

function increment(i,c,N)
  if i == 1
    disp(c);
    while c(i) < N(i)
      c(i) = c(i) + 1;
      disp(c);
    endwhile
  elseif i > 1
    while c(i) <= N(i)
      increment(i-1,c,N);
      c(i) = c(i) + 1;
    endwhile
  endif
endfunction

increment(size(R)(2),d,R);