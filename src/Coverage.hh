//
// varViewer
// Date: 2023-08-24
// Author: Huang Fei <huangfei@genomics.cn>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#pragma once

#include <iostream>

namespace varViewer {

//static const size_t COVS_ALIGN = 7;
/* enum covs_idx {
*    SNP_VR = 0,  // SNP variant reads
*    SNP_RR = 1,  // SNP reference reads
*    INS_VR = 2,  // Insert variant reads
*    INS_RR = 3,  // Insert reference reads
*    DEL_VR = 4, //  Deletion variant reads
*    DEL_RR = 5,  // Deletion reference reads
*    INSDEL_VR = 6 // Insdel variant reads
* };
*/

static const size_t COVS_ALIGN = 2;

enum covs_idx {
    VAR_VR =0, // SNP/InDel variant reads
    VAR_RR = 1 // SNP/InDel reference reads
};

class Coverage {
  private:
    int (*_covs)[COVS_ALIGN];
    int _len;

    Coverage();

  public:
    Coverage(size_t len);
    ~Coverage();
    void reset();
    void add_coverage(int pos, covs_idx idx);
    void add_coverage(int pos, covs_idx idx, int add_cov);
    void add_coverage_range(int pos, covs_idx idx, int len);
    void add_coverage_range(int pos, covs_idx idx, int len, int add_cov);
    void set_coverage(int pos, covs_idx idx, int cov);
    int get_coverage(int pos, covs_idx idx);
    int get_pos_coverage(int pos);
    int get_max_coverage();
    int get_len(){ return this->_len;}
    int get_min_coverage(int start, int end, covs_idx idx);
};
std::ostream& operator<<(std::ostream& os, Coverage& cover);
}