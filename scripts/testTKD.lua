-- creates a ugx file
-- author: Martin Scherer

-- parameter fuer einheits tkd
a = 10
w = 3*a
h = math.sqrt(6)*a 
d_lipid = 5

rows = 4 
cols = 2 
high = 2 

--rows = 4
--cols = 2
--high = 2

--rows = 1
--cols = 1
--high = 1

date = os.date("%d-%m-%y__%H-%M-%S")
os.execute("mkdir -p /tmp/tkd")
filename = "/tmp/tkd/tkd__" .. date .. ".ugx"

-- calling testing method
TestTKDGenerator(filename, h, a, w, d_lipid, rows, cols, high)

print("Grid written to " .. filename)