function c=hex2rgb(s)

c=[hex2dec(s(1:2)) hex2dec(s(3:4)) hex2dec(s(5:6))]/255;
