
row_perm = "0 8 6 1 15 10 14 7 2 3 12 13 4 11 9 5"
col_perm = "15 9 2 0 5 7 11 14 3 12 1 6 10 8 13 4"

slp = """y00 = x03 + x09 + x15
y01 = x03 + x13 + y00
y02 = x02 + x11 + x12
y03 = x04 + x10 + x14
y04 = x00 + x06 + x08
y05 = x02 + x07 + y02
y06 = x02 + x07 + x12
y07 = x00 + x05 + y04
y08 = x02 + x07 + x11
y09 = x00 + x05 + x08
y10 = x03 + x09 + x13
y11 = x01 + x04 + y03
y12 = x03 + x09 + y01
y13 = x00 + x05 + x06
y14 = x01 + x04 + x14
y15 = x01 + x04 + x10"""

col_dict = {"x" + str(k).zfill(2): "x" + str(v).zfill(2) for k, v in enumerate(col_perm.split(" "))}
row_dict = {"y" + str(k).zfill(2): "y" + str(v).zfill(2) for k, v in enumerate(row_perm.split(" "))}

for line in slp.split("\n"):
    lhs, rhs = line.strip().split(' = ')
    
    if 'x' in lhs:
        lhs = col_dict[lhs]
    elif 'y' in lhs:
        lhs = row_dict[lhs]
    
    rhs_new = []

    if '+' not in rhs:
        if 'x' in rhs:
            rhs_new.append(col_dict[rhs])
        elif 'y' in rhs:
            rhs_new.append(row_dict[rhs])
    else:
        for i in rhs.split(' + '):
            if 'x' in i:
                rhs_new.append(col_dict[i])
            elif 'y' in i:
                rhs_new.append(row_dict[i])
            else:
                rhs_new.append(i)
    
    if len(rhs_new) == 1:
        print(lhs + " = " + rhs_new[0])
    else:
        print(lhs + " = " + ' + '.join(rhs_new))
