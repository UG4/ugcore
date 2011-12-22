-- Creates a single tkd
-- author: Martin Scherer

-- einheits tkd
a = 10
w = 3*a
h = math.sqrt(6)*a 
d_lipid = 8

rows = 1
cols = 1
high = 1

--rows= 3
--cols =4
--high=2

date = os.date("%d-%m-%y__%H-%M-%S")
filename = "tkd__" .. date .. ".ugx"

-- calling testing method
TestTKDGenerator(filename, h, a, w, d_lipid, rows, cols, high)

print ("Grid written to " .. filename)