--
-- -----------------------------------------------------------------
-- COMPANY : Ruhr University Bochum
-- AUTHOR  : Amir Moradi amir.moradi@rub.de Aein Rezaei Shahmirzadi aein.rezaeishahmirzadi@rub.de
-- DOCUMENT: "Second-Order SCA Security with almost no Fresh Randomness" TCHES 2021, Issue 3
-- -----------------------------------------------------------------
--
-- Copyright c 2021, Amir Moradi, Aein Rezaei Shahmirzadi
--
-- All rights reserved.
--
-- THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
-- ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
-- WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
-- DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTERS BE LIABLE FOR ANY
-- DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
-- INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
-- LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION HOWEVER CAUSED AND
-- ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
-- INCLUDING NEGLIGENCE OR OTHERWISE ARISING IN ANY WAY OUT OF THE USE OF THIS
-- SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
--
-- Please see LICENSE and README for license and further instructions.
--

LIBRARY ieee;
USE ieee.std_logic_1164.ALL;
USE ieee.numeric_std.ALL;
USE ieee.math_real.ALL; 

ENTITY PRINCE_TB IS
END PRINCE_TB;
 
ARCHITECTURE behavior OF PRINCE_TB IS 

   --Inputs
   signal input_s1 :	std_logic_vector(63 downto 0) := (others => '0');
   signal input_s2 :	std_logic_vector(63 downto 0) := (others => '0');
   signal input_s3 :	std_logic_vector(63 downto 0) := (others => '0');
   signal Key1 :			std_logic_vector(127 downto 0) := (others => '0');
   signal Key2 :			std_logic_vector(127 downto 0) := (others => '0');
   signal Key3 :			std_logic_vector(127 downto 0) := (others => '0');
   signal random :			std_logic_vector(127 downto 0) := (others => '0');
	signal enc_dec :		std_logic := '0';
   signal clk :			std_logic := '0';
   signal reset :			std_logic := '0';

 	--Outputs
   signal output_s1 :	std_logic_vector(63 downto 0);
   signal output_s2 :	std_logic_vector(63 downto 0);
   signal output_s3 :	std_logic_vector(63 downto 0);
   signal done :			std_logic;

	signal input  :	std_logic_vector(63 downto 0) := (others => '0');
	signal output :	std_logic_vector(63 downto 0);
   
	-- Clock period definitions
   constant clk_period : time := 10 ns;
 
	type INT_ARRAY  is array (integer range <>) of integer;
	type REAL_ARRAY is array (integer range <>) of real;
	type BYTE_ARRAY is array (integer range <>) of std_logic_vector(7 downto 0);

	signal r: INT_ARRAY (15 downto 0);
	signal m: BYTE_ARRAY(15 downto 0);
 
BEGIN
 
  	maskgen: process
		 variable seed1, seed2: positive;        -- seed values for random generator
		 variable rand: REAL_ARRAY(15 downto 0); -- random real-number value in range 0 to 1.0  
		 variable range_of_rand : real := 256.0; -- the range of random values created will be 0 to +1000.
	begin
		 
		FOR i in 0 to 15 loop
			uniform(seed1, seed2, rand(i));   -- generate random number
			r(i) <= integer(rand(i)*range_of_rand);  -- rescale to 0..1000, convert integer part 
			m(i) <= std_logic_vector(to_unsigned(r(i), m(i)'length));
		end loop;
		
		wait for clk_period;
	end process;  

	---------

--	input_s2 <= m(0)  & m(1)  & m(2)  & m(3)  & m(4)  & m(5)  & m(6)  & m(7);
--	input_s3 <= m(8)  & m(9)  & m(10) & m(11) & m(12) & m(13) & m(14) & m(15);
 	input_s2 <= (others => '0');
	input_s3 <= (others => '0');
	---------
	
	output <= output_s1 xor output_s2 xor output_s3;
  
	-- Instantiate the Unit Under Test (UUT)
   uut: entity work.PRINCE 
	PORT MAP
	(
		clk => clk,
           reset => reset,
           input_s1 => input_s1,
           input_s2 => input_s2,
           input_s3 => input_s3,
           r => random,
           output_s1 => output_s1,
           output_s2 => output_s2,
           output_s3 => output_s3,
			  Key1 => Key1,
			  Key2 => Key2,
			  Key3 => Key3,
			  enc_dec => enc_dec,
			  done => done
	);

	input_s1 <= input ;--xor input_s2 xor input_s3;

   -- Clock process definitions
   clk_process :process
   begin
		clk <= '1';
		wait for clk_period/2;
		clk <= '0';
		wait for clk_period/2;
   end process;
 

   -- Stimulus process
   stim_proc: process
   begin		
      -- hold reset state for 100 ns.
--      wait for 100 ns;	
		wait for clk_period/2;

			input 	<= x"0123456789ABCDEF"; 
			Key1 		<= x"FEDCBA9876543210FEDCBA9876543210";
			Key2		<= (others => '0');
			Key3		<= (others => '0');
			enc_dec 	<= '0';
			
			wait for clk_period;
			
			reset <= '1';
			
			wait for 7*clk_period;
					
			reset <= '0';
			
			wait until done='1';
			
			wait for 0.5*clk_period;
			
			if (output = x"D9830DF8619840CC") then
				report "---------- Passed ----------";
			else
				report "---------- Failed ----------";
			end if;	
			
			wait for 10*clk_period;

			input 	<= x"D9830DF8619840CC";
			Key1 		<= x"FEDCBA9876543210FEDCBA9876543210";
			enc_dec 	<= '1';
			
			wait for clk_period;
			
			reset <= '1';
			
			wait for 8*clk_period;
					
			reset <= '0';

			wait until done='1';
			
			wait for 0.5*clk_period;
			
			if (output = x"0123456789ABCDEF") then
				report "---------- Passed ----------";
			else
				report "---------- Failed ----------";
			end if;	
			
			wait for 10*clk_period;

      wait;
   end process;

END;