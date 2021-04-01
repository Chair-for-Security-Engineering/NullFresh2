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

-- IMPORTS
----------------------------------------------------------------------------------
LIBRARY IEEE;
USE IEEE.STD_LOGIC_1164.ALL;



-- ENTITY
----------------------------------------------------------------------------------
ENTITY TB_Skinny IS
END TB_Skinny;



-- ARCHITECTURE : BEHAVIOR
----------------------------------------------------------------------------------
ARCHITECTURE Behavioral OF TB_Skinny IS

	-- CONSTANTS ------------------------------------------------------------------
	CONSTANT N : INTEGER := 64;
	CONSTANT T : INTEGER := 64;
	CONSTANT W : INTEGER := 4;
	-------------------------------------------------------------------------------

	-- TEST VECTORS ---------------------------------------------------------------
	SIGNAL TV_PT 	: STD_LOGIC_VECTOR((N - 1) DOWNTO 0);
	SIGNAL TV_CT 	: STD_LOGIC_VECTOR((N - 1) DOWNTO 0);
	SIGNAL TV_KT 	: STD_LOGIC_VECTOR((T - 1) DOWNTO 0);
	SIGNAL M1 		: STD_LOGIC_VECTOR(63      DOWNTO 0);
	SIGNAL M2 		: STD_LOGIC_VECTOR(63      DOWNTO 0);
	SIGNAL M3 		: STD_LOGIC_VECTOR(63      DOWNTO 0);
	SIGNAL M4 		: STD_LOGIC_VECTOR(63      DOWNTO 0);
	
	-------------------------------------------------------------------------------

	-- INPUTS ---------------------------------------------------------------------
   SIGNAL CLK 			: STD_LOGIC := '0';
   SIGNAL RESET 		: STD_LOGIC := '0';
   SIGNAL FRESH      : STD_LOGIC_VECTOR(127     DOWNTO 0) := (others => '0');
	SIGNAL KEY1 		: STD_LOGIC_VECTOR((T - 1) DOWNTO 0);
   SIGNAL KEY2 		: STD_LOGIC_VECTOR((T - 1) DOWNTO 0);
   SIGNAL KEY3 		: STD_LOGIC_VECTOR((T - 1) DOWNTO 0);
   SIGNAL PLAINTEXT1	: STD_LOGIC_VECTOR((N - 1) DOWNTO 0);
   SIGNAL PLAINTEXT2	: STD_LOGIC_VECTOR((N - 1) DOWNTO 0);
   SIGNAL PLAINTEXT3	: STD_LOGIC_VECTOR((N - 1) DOWNTO 0);
	-------------------------------------------------------------------------------

	-- OUTPUTS --------------------------------------------------------------------
	SIGNAL DONE			: STD_LOGIC;
   SIGNAL CIPHERTEXT : STD_LOGIC_VECTOR((N - 1) DOWNTO 0);
   SIGNAL CIPHERTEXT1: STD_LOGIC_VECTOR((N - 1) DOWNTO 0);
   SIGNAL CIPHERTEXT2: STD_LOGIC_VECTOR((N - 1) DOWNTO 0);
   SIGNAL CIPHERTEXT3: STD_LOGIC_VECTOR((N - 1) DOWNTO 0);
	-------------------------------------------------------------------------------

   -- CLOCK PERIOD DEFINITIONS ---------------------------------------------------
   CONSTANT CLK_PERIOD : TIME := 10 NS;
	-------------------------------------------------------------------------------

BEGIN

	-- INSTANTIATE UNIT UNDER TEST (UUT) ------------------------------------------
   UUT : ENTITY work.Skinny
	PORT MAP (
		CLK 			=> CLK,
		RESET 		=> RESET,
		FRESH			=> FRESH,
		DONE			=> DONE,
		KEY1 			=> KEY1,
		KEY2 			=> KEY2,
		KEY3 			=> KEY3,
		PLAINTEXT1 	=> PLAINTEXT1,
		PLAINTEXT2 	=> PLAINTEXT2,
		PLAINTEXT3 	=> PLAINTEXT3,
		CIPHERTEXT1	=> CIPHERTEXT1,
		CIPHERTEXT2	=> CIPHERTEXT2,
		CIPHERTEXT3	=> CIPHERTEXT3);
	-------------------------------------------------------------------------------

	CIPHERTEXT <= CIPHERTEXT1 XOR CIPHERTEXT2 XOR CIPHERTEXT3;

   -- CLOCK PROCESS --------------------------------------------------------------
   CLK_PROCESS : PROCESS
	BEGIN
		CLK <= '0'; WAIT FOR CLK_PERIOD/2;
		CLK <= '1'; WAIT FOR CLK_PERIOD/2;
   END PROCESS;
	-------------------------------------------------------------------------------

   -- STIMULUS PROCESS -----------------------------------------------------------
   STIM_PROCESS : PROCESS
   BEGIN

		----------------------------------------------------------------------------
		TV_PT <= X"DC6A4944FF7DE341";
		TV_CT <= X"6DE8CEA4CAF93063";
		TV_KT	<= X"E8619935055EE7BF";
		
		M1		<= X"DB502D1AE642C60B";
		M2		<= X"7F292E3E59BFCA6F";
		M3		<= X"209542CF81B9F2D6";
		M4		<= X"77854D5D8A150564";
				
		----------------------------------------------------------------------------

		WAIT FOR CLK_PERIOD;

		----------------------------------------------------------------------------
		PLAINTEXT1	<= TV_PT XOR M1 XOR M2;
		PLAINTEXT2	<= M1;
		PLAINTEXT3	<= M2;
		
		KEY1			<= TV_KT XOR M3 XOR M4;
		KEY2			<= M3;
		KEY3			<= M4;
		----------------------------------------------------------------------------

		----------------------------------------------------------------------------
		RESET <= '1';
			WAIT FOR CLK_PERIOD*4;
		RESET <= '0';
		----------------------------------------------------------------------------

     	WAIT UNTIL DONE = '1';

		----------------------------------------------------------------------------
	   WAIT FOR CLK_PERIOD/2;
		----------------------------------------------------------------------------

  	 	----------------------------------------------------------------------------
     	IF (CIPHERTEXT = TV_CT) THEN
   		ASSERT FALSE REPORT "---------- PASSED ----------" SEVERITY WARNING;
     	ELSE
         ASSERT FALSE REPORT "---------- FAILED ----------" SEVERITY WARNING;
    	END IF;
		----------------------------------------------------------------------------
		
		WAIT FOR CLK_PERIOD*100;

   END PROCESS;
	-------------------------------------------------------------------------------

END;
