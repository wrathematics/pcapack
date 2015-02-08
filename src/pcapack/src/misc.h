/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2015, Schmidt

#ifndef __PCAPACK_MISC_H__
#define __PCAPACK_MISC_H__


#define MIN(m,n) m<n?m:n
#define MAX(m,n) m<n?n:m

#define likely(x)   __builtin_expect((x),1)
#define unlikely(x) __builtin_expect((x),0)


#endif


