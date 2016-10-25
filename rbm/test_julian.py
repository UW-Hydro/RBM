import leapyear
import numpy as np

test_date = [1900, 2000, 2001, 2004]
answer = [0,1,0,1]

for i in range(len(test_date)):
    determine = leapyear.leapyear(test_date[i]) == answer[i]
    if determine != True:
        print('Wrong!')
