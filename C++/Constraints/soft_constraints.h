// Copyright 2019 Alexander Liniger

// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at

//     http://www.apache.org/licenses/LICENSE-2.0

// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

#ifndef MPCC_SOFTCONSTRAINTS_H
#define MPCC_SOFTCONSTRAINTS_H


#include "config.h"
#include "types.h"
namespace mpcc {
class SoftConstraints {

private:
    Q_MPC Zx;
    R_MPC Zu;

    q_MPC zx;
    r_MPC zu;
};
}

#endif //MPCC_SOFTCONSTRAINTS_H
