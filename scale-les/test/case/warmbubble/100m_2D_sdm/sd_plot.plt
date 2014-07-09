reset
st = 0
to = 1790
offset = 10
do for[i = st : to : offset]{
    set xrange [0:1000]
    set yrange [0:400]
    set zrange [0:4000]
    input0 = sprintf("SD_output_ASCII_0000000%04d.000.pe000000", i)
    input1 = sprintf("SD_output_ASCII_0000000%04d.000.pe000001", i)
    input2 = sprintf("SD_output_ASCII_0000000%04d.000.pe000002", i)
    input3 = sprintf("SD_output_ASCII_0000000%04d.000.pe000003", i)
#    plot input0 using 2:3, input1 using 2:3, input2 using 2:3, input3 using 2:3
#    plot input0 using 1:3, input1 using 1:3, input2 using 1:3, input3 using 1:3
    plot input0 using 1:( $3>0 ? $2 : 1/0), input1 using 1:( $3>0 ? $2 : 1/0), input2 using 1:( $3>0 ? $2 : 1/0), input3 using 1:( $3>0 ? $2 : 1/0)
#    splot input0, input1, input2, input3
    pause 0.5
}
reset
