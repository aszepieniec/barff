F = FiniteField(4096)
z = F.gen()

def F2hex( elm ):
    coeffs = elm.polynomial().coefficients(sparse=False)
    integer = sum(2^i * ZZ(coeffs[i]) for i in range(0,len(coeffs)))
    return hex(integer)

def int2F( integer ):
    bits = bin(integer)[2:]
    bits = [bits[len(bits)-1-i] for i in range(0, len(bits))]
    return sum(F(bits[i]) * z^i for i in range(0, len(bits)) if bits[i] == '1')

def dlog( elm ):
    if elm == 0:
        return 4095
    for i in range(0, 4095):
        if z^i == elm:
            return i

def antilog( integer ):
    return z^integer

print "#ifndef GF4096_TABLES_H"
print "#define GF4096_TABLES_H"
print ""
print "int gf4096_antilogs[4095] = {"
for i in range(0, 4095):
    print "0x%s" % F2hex(z^i),
    if i == 4094:
        print "};"
    elif i % 8 == 7:
        print ","
    else:
        print ",",
print ""
print "int gf4096_dlogs[4096] = {"
for i in range(0, 4096):
    print dlog(int2F(i)),
    if i == 4095:
        print "};"
    elif i % 8 == 7:
        print ","
    else:
        print ",",
print ""
print "#endif"

