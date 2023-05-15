import math
a=float(input())
b=float(input())
op=input("what do you want?")
if op!=None:
    if op=="+":
         result=a+b
    elif op=="-":
         result=a-b
    elif op=="*":
        result=a*b
    elif op=="/":
        if b==0:
            result="error"
        else:
            result=a/b
    elif op=="sin":
        result=("sin a=",math.sin(a),"sin b=",math.sin(b))
    elif op=="tan":
        result=("tan a=",math.tan(a),"tan b=", math.tan(b))
    elif op=="tan":
        result=("tan a=",math.tan(a),"tan b=", math.tan(b))
    print(result)
else:
    print("invalid")


