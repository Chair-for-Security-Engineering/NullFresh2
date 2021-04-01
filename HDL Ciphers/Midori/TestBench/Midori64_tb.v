`timescale 1ns / 1ps

/*
* -----------------------------------------------------------------
* COMPANY : Ruhr University Bochum
* AUTHOR  : Amir Moradi amir.moradi@rub.de Aein Rezaei Shahmirzadi aein.rezaeishahmirzadi@rub.de
* DOCUMENT: "Second-Order SCA Security with almost no Fresh Randomness" TCHES 2021, Issue 3
* -----------------------------------------------------------------
*
* Copyright c 2021, Amir Moradi, Aein Rezaei Shahmirzadi
*
* All rights reserved.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTERS BE LIABLE FOR ANY
* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* INCLUDING NEGLIGENCE OR OTHERWISE ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
* Please see LICENSE and README for license and further instructions.
*/

module Midori64_tb;

	// Inputs
	reg clk;
	reg reset;
	reg [63:0] input1;
	reg [63:0] input2;
	reg [63:0] input3;
	reg [127:0] Key1;
	reg [127:0] Key2;
	reg [127:0] Key3;
	reg [127:0] r;
	reg enc_dec;


	// Outputs
	wire [63:0] output1;
	wire [63:0] output2;
	wire [63:0] output3;
	
	
	wire [63:0] Out;
	wire [63:0] In;
	wire done;

	assign In = input1 ^ input2 ^ input3;
	assign Out = output1 ^ output2 ^ output3;
	// Instantiate the Unit Under Test (UUT)
	Midori64 uut (
		.clk(clk), 
		.reset(reset), 
		.input1(input1), 
		.input2(input2), 
		.input3(input3), 
		.output1(output1), 
		.output2(output2), 
		.output3(output3), 
		.Key1(Key1), 
		.Key2(Key2), 
		.Key3(Key3), 
		.r(r), 
		.enc_dec(enc_dec), 
		.done(done)
	);

	initial begin
		// Initialize Inputs
		clk = 0;
		reset = 1;
		input1 = 0;
		input2 = 0;
		input3 = 0;
		Key1 = 128'h687ded3b3c85b3f35b1009863e2a8cbf;
		Key2 = 0;
		Key3 = 0;
		enc_dec = 0;
		r = 0;

		// Wait 100 ns for global reset to finish
		#10;
		input1 = 64'h42c20fd3b586879e;
		#10;
		input1 = 64'h0;
		#10;
		input1 = 64'h42c20fd3b586879e;
		#10;
		reset = 0;
		// Add stimulus here
		@(posedge done) begin
			#5;
			if(Out == 64'h36f32dcf124ab057) begin
				$write("------------------PASS---------------\n");
			end
			else begin
				$write("\------------------FAIL---------------\n");
				$write("%x\n%x\n%x\n", In, Out, 64'h36f32dcf124ab057);
			end
			#10;
			if(Out == 64'h66bcdc6270d901cd) begin
				$write("------------------PASS---------------\n");
			end
			else begin
				$write("\------------------FAIL---------------\n");
				$write("%x\n%x\n%x\n", In, Out, 64'h66bcdc6270d901cd);
			end
			#10;
			if(Out == 64'h36f32dcf124ab057) begin
				$write("------------------PASS---------------\n");
			end
			else begin
				$write("\------------------FAIL---------------\n");
				$write("%x\n%x\n%x\n", In, Out, 64'h36f32dcf124ab057);
			end
			#10;
			if(Out == 64'h66bcdc6270d901cd) begin
				$write("------------------PASS---------------\n");
			end
			else begin
				$write("\------------------FAIL---------------\n");
				$write("%x\n%x\n%x\n", In, Out, 64'h66bcdc6270d901cd);
			end
		end
		#10;
		$stop;

	end
      always #5 clk=~clk;
endmodule

