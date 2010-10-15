print("Look at my new cake!")

c = Cake()
p = c:take_pieces(7)

print("If i take " .. p:size() .. " pieces from it, there are still " .. c:pieces_left() .. " left.")
print("now I'll put them back.")

c:add_pieces(p)

print("And now there are " .. c:pieces_left() .. " pieces.")
print("I can't understand this.")

print("")
print("Creating base and derived.")
print("The function 'PrintFunction(Base& b)' is only callable with Base and Derived. We test that now ... \n")

-- base = Base()
derived = Derived()

print("Calling PrintFunction with Derived:")
PrintFunction(derived)
print("\nCalling PrintFunction with Base:")
-- PrintFunction(base)
print("\nCalling PrintFunction with Cake:")
PrintFunction(c)
print("\nCalling PrintFunction with 7:")
PrintFunction(7)

print("done")
print("Great!")
