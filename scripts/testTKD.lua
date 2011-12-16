-- Creates a single tkd
-- author: Martin Scherer

-- einheits tkd
a = 10
w = 3*a
h = math.sqrt(6)*a
d_lipid = 9

date = os.date("%d-%m-%y__%H-%M-%S")
filename = "/tmp/tkd__" .. date .. ".ugx"

-- calling testing method
TestTKDGenerator(filename, h, a, w, d_lipid)

print ("Grid written to " .. filename)