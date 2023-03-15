open("data/VenturaAccelerometer/test_mesh_01.msh") do f
 
    # line_number
    line = 0  
   
    # read till end of file
    while ! eof(f) 
   
       # read a new / next line for every iteration          
       s = readline(f)         
       line += 1
       println("$line . $s")
       if ("$line . $s") == ("\$Nodes")
        print("HEYYYY")
       end
    end
   
  end
