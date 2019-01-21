/****************************************************************/
/*                                                              */
/*  Copyright (C) 2018, M. Andelkovic, L. Covaci, A. Ferreira,  */
/*                    S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                              */
/****************************************************************/
#include "Generic.hpp"
#include "ComplexTraits.hpp"
#include "myHDF5.hpp"
#include <algorithm>
#include "queue.hpp"

std::string num2str3(int dir_num){
  std::string dir;
  
  switch(dir_num){
  case 0:
    dir = "xxx"; break;
  case 1:
    dir = "xxy"; break;
  case 2:
    dir = "xxz"; break;
  case 3:
    dir = "xyx"; break;
  case 4:
    dir = "xyy"; break;
  case 5:
    dir = "xyz"; break;
  case 6:
    dir = "xzx"; break;
  case 7:
    dir = "xzy"; break;
  case 8:
    dir = "xzz"; break;
  case 9:
    dir = "yxx"; break;
  case 10:
    dir = "yxy"; break;
  case 11:
    dir = "yxz"; break;
  case 12:
    dir = "yyx"; break;
  case 13:
    dir = "yyy"; break;
  case 14:
    dir = "yyz"; break;
  case 15:
    dir = "yzx"; break;
  case 16:
    dir = "yzy"; break;
  case 17:
    dir = "yzz"; break;
  case 18:
    dir = "zxx"; break;
  case 19:
    dir = "zxy"; break;
  case 20:
    dir = "zxz"; break;
  case 21:
    dir = "zyx"; break;
  case 22:
    dir = "zyy"; break;
  case 23:
    dir = "zyz"; break;
  case 24:
    dir = "zzx"; break;
  case 25:
    dir = "zzy"; break;
  case 26:
    dir = "zzz"; break;
  default:
    std::cout << "Invalid direction in num2str_dir3.\n"; exit(1);
  }
  return dir;
}

std::string num2str2(int dir_num) {
  std::string dir;
  
  switch(dir_num){
  case 0:
    dir = "xx"; break;
  case 1:
    dir = "yy"; break;
  case 2:
    dir = "zz"; break;
  case 3:
    dir = "xy"; break;
  case 4:
    dir = "xz"; break;
  case 5:
    dir = "yx"; break;
  case 6:
    dir = "yx"; break;
  case 7:
    dir = "zx"; break;
  case 8:
    dir = "zy"; break;
  default:
    std::cout << "Invalid direction for the optical conductivity.\n"; exit(1);
  }
  return dir;
}

