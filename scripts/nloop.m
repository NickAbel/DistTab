R = [2,2,2];
d = [1,1,1];

function increment(i,c,N)
  if i == 1
    disp(flip(c));
    while c(i) < N(i)
      c(i) = c(i) + 1;
      disp(flip(c));
    endwhile
  elseif i > 1
    while c(i) <= N(i)
      increment(i-1,c,N);
      c(i) = c(i) + 1;
    endwhile
  endif
endfunction

increment(size(R)(2),d,R);