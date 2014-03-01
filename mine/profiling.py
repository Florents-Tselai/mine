import pandas as pd

d = {
    'Linear':"$y=2^x$",
    'Parabolic':"$y=2*x+1$",
    'Cubic':"$y=3*x^2+2*x+1$",
    'Exponential':"$y=10^(10x)-1$",
    'Linear/Periodic':"y=sin(10pix)+x",
    'Sinusodial (Fourier Frequency)':"y=sin(16pix)",
    'Sinusodial (non-Fourier Frequency':"y=sin(13pix)",
    'Sinusodial (Varying Frequency':"y=sin(7pix(1+x)",
    'Categorical':"64 points chosen from: {}",
    'Random':"random number generator"
}

df = pd.DataFrame(columns=['Relationship Name', 'Description', 100, 1000, 10000, 100000])
df['Relationship Name'] = d.keys()
#df.set_index('Relationship Name', inplace=True)
df['Description'] =  df['Relationship Name'].map(lambda x: d[x])
print df.index
output = df.to_latex()
output = output.replace("\$", "$").replace("\\textasciicircumx", "").replace("\\textasciicircum", "")
print output
text_file = open("/home/florents/workspace/mine/doc/profiling_1.tex", "w")
text_file.write(output)
text_file.close()
