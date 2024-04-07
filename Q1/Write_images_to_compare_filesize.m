
imwrite(uint8(rec),"rec.jpg")
imwrite(imagefirst,"first.jpg")

first = dir('first.jpg');         
fs = first.bytes

rec = dir('rec.jpg');         
rs = rec.bytes

compression=100-(rs/fs)*100