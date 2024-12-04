/***********************************************************************[version.hpp]
Copyright(c) 2020, Muhammad Osama - Anton Wijs,
Technische Universiteit Eindhoven (TU/e).

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
**********************************************************************************/

namespace ParaFROST {

	const char* version();
	const char* signature();
	const char* compiler();
	const char* compilemode();
	const char* osystem();
	const char* date();

}
#define VERSION "3.4.1"
#define OSYSTEM "linux asus-f16 5.15.167.4-microsoft-standard-wsl2 x86_64"
#define DATE "Wed Dec  4 14:33:41 EST 2024"
#define COMPILER "g++ (Ubuntu 11.4.0-1ubuntu1~22.04) 11.4.0 + nvcc 12.6.85"
