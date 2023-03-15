open("data/VenturaAccelerometer/test_mesh_01.msh") do f
 
    # line_number
    line = 0  

NON = 0  
NOE = 0
fnode = 0
enode = 0
felement = 0
eelement = 0

    # read till end of file
    while ! eof(f) 
   
       # read a new / next line for every iteration          
       s = readline(f)         
       line += 1
       println("$line . $s")
        if "$s" == "\$Nodes"
            NON = parse(Int, "$line") + 1
            fnode = NON + 1
            print(NON)
            print("-----------------\n")
            print(fnode)
            print("-----------------\n")
        elseif "$s" == "\$EndNodes"
            enode = parse(Int, "$line") -1
            print(enode)
            print("-----------------\n")
        elseif "$s" == "\$Elements"
            NOE = parse(Int, "$line") + 1
            felement = NOE + 1
            print(NOE)
            print("-----------------\n")
            print(felement)
            print("-----------------\n")
        elseif "$s" == "\$EndElements"
            eelements = parse(Int, "$line") - 1
            print("-----------------\n")
            break
        end
    end

    
end
