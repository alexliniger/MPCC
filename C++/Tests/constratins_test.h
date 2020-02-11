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

#ifndef MPCC_CONSTRATINS_TEST_H
#define MPCC_CONSTRATINS_TEST_H


#include "Spline/arc_length_spline.h"
#include "Constraints/constraints.h"
namespace mpcc{
void genRoundTrack(ArcLengthSpline& track);
int testAlphaConstraint(const PathToJson &path);
int testTireForceConstraint(const PathToJson &path);
int testTrackConstraint(const PathToJson &path);
}
#endif //MPCC_CONSTRATINS_TEST_H
