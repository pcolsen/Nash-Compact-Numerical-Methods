#!/usr/bin/python3
# Code to Measure time taken by program to execute. 
import time
import math

n = 2345678

# store starting time 
begin = time.time() 

# program body starts
y = 0.0
for i in range(n):
   x = math.exp(math.sin(math.cos(i)))
   y = y+x

print("Sum of outputs=",y)
# program body ends 

time.sleep(1) 
# store end time 
end = time.time() 

# total time taken 
print(f"Total runtime of the program is {end - begin}") 
