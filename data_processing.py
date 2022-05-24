import pandas as pd
import matplotlib.pyplot as plt

def read_table_pandas_x(filename):
    data_x = pd.read_csv(filename, skip_blank_lines=True, delim_whitespace=True, header=None, usecols=[2])
    return data_x

def read_table_pandas_y(filename):
    data_y = pd.read_csv(filename, skip_blank_lines=True, delim_whitespace=True, header=None, usecols=[0])
    return data_y

data_x = read_table_pandas_x('./Cases/Bayesian/cerbereoutput/probe_EQ.out')
data_y = read_table_pandas_y('./Cases/Bayesian/cerbereoutput/probe_EQ.out')
print(data_y)
plt.plot(data_y, data_x)
plt.semilogx()
plt.show()
