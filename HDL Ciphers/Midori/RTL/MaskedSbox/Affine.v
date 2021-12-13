module Affine(
    input [3:0] x1,
    input [3:0] x2,
    input [3:0] x3,
    output [3:0] y1,
    output [3:0] y2,
    output [3:0] y3
    );



	parameter num = 1;
	
	wire [3:0] notx1;
	assign notx1 = ~x1;
	
	generate
		if(num == 1) begin
			assign y1 = {notx1[0] ^ x1[2],	x1[3],	x1[0] ^ x1[3],	notx1[1]};
			assign y2 = {x2[0] ^ x2[2], 	x2[3], 	x2[0] ^ x2[3], 	x2[1]};
			assign y3 = {x3[0] ^ x3[2], 	x3[3], 	x3[0] ^ x3[3], 	x3[1]};

		end
   endgenerate

endmodule
