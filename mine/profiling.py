import pandas as pd

d = {
    1:"$f(x)=2^x$",
    2:"$f(x)=2*x+1$",
    3:"$f(x)=3*x^2+2*x+1$"
}

df = pd.DataFrame(index=d.keys(), columns=['function', "100", "200", "300"])
#print df.index
df['function'] = map(lambda x: d[x], df.index)

df = df.T
output = df.to_latex()
output = output.replace("\$", "$").replace("\\textasciicircumx", "").replace("\\textasciicircum", "")
print output
text_file = open("/home/florents/workspace/mine/doc/profiling_1.tex", "w")
text_file.write(output)
text_file.close()
