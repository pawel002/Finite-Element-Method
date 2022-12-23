import matplotlib.pyplot as plt

plt.grid(visible=True)

with open('.\\saveplot\\plot.txt') as f:
    lines = f.readlines()

X = []
Y = []

for line in lines:
    line.replace('\n', '')
    x = line.split(' ')
    a = float(x[0])
    b = float(x[1])

    X.append(a)
    Y.append(b)

plt.plot(X, Y, "blue")
plt.show()
