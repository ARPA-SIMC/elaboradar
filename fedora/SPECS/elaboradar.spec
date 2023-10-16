# Note: define srcarchivename in CI build only.
%{!?srcarchivename: %global srcarchivename %{name}-%{version}-%{release}}

Summary:	Library and tools to handle weather radar images and data
Name: 		elaboradar
Version: 	0.23
Release: 	1
License: 	GPL
Group: 		Applications/Meteo
URL:            https://github.com/arpa-simc/%{name}
Source:         https://github.com/arpa-simc/%{name}/archive/v%{version}-%{release}.tar.gz#/%{srcarchivename}.tar.gz

BuildRoot:      %{_tmppath}/%{name}-%{version}-%{release}-buildroot
Packager: 	Daniele Branchini <dbranchini@arpae.it>
BuildRequires:	libtool, gcc-c++, hdf5-devel, eigen3-devel, radarlib-devel, log4c-devel, gsl-devel, tclap, gdal-devel, proj-devel, doxygen
Requires:       hdf5

%description
Library and to handle weather radar images and data

%package devel
Requires: elaboradar = %{version}
Group: Libraries/Meteo
Summary: Development for radarelab library

%description devel
Development for elaboradar library

%package doc
Summary: elaboradar documentation
Group: Libraries/Meteo

%description doc
elaboradar library documentation

%package tools
Requires: elaboradar = %{version}
Group: Applications/Meteo
Summary: Development for radarelab library

%description tools
Tools for elaboradar library

%prep
%setup -q -n %{srcarchivename}
sh autogen.sh

%build

%configure
make
make check

%install
[ "%{buildroot}" != / ] && rm -rf %{buildroot}
%makeinstall

%clean
[ "%{buildroot}" != / ] && rm -rf %{buildroot}

%files
%defattr(-,root,root,-)
%{_libdir}/libradarelab.so.0*

%files devel
%defattr(-,root,root,-)
%{_includedir}/radarelab/*
%{_libdir}/libradarelab.a
%exclude %{_libdir}/libradarelab.la
%{_libdir}/libradarelab.so
%{_libdir}/pkgconfig/radarelab.pc

%files doc
%defattr(-,root,root,-)
%doc %{_docdir}/%{name}

%files tools
%{_bindir}/AddCleanerQuantities
%{_bindir}/classificatore
%{_bindir}/elaboradar
%{_bindir}/RunCleanID
%{_bindir}/stat_CleanID
%dir %{_libexecdir}/%{name}
%{_libexecdir}/%{name}/vecchioripulisco
%{_datadir}/%{name}/*

%changelog
* Mon Oct 16 2023 Daniele Branchini <dbranchini@arpae.it> - 0.23-1
- Corretto problema buffer overflow su stat_CleanID

* Fri Oct 13 2023 Daniele Branchini <dbranchini@arpae.it> - 0.22-1
- Rimosso il nome radar dagli argomenti passati in input da riga di comando a RunCleanID e stat_CleanID
- Aggiunto il numero di elevazioni da stampare tra gli argomenti passati a stat_CleanID
- Aggiunta la stampa del numero totale di bin per ogni raggio per ogni elevazione (beam_size) in stat_CleanID

* Tue Jun 13 2023 Daniele Branchini <dbranchini@arpae.it> - 0.21-1
- Reduced test verbosity

* Mon Jun 12 2023 Daniele Branchini <dbranchini@arpae.it> - 0.20-1
- Aggiunta mappe statiche corrette su 2 settori di occlusione valide fino al 30/5/2023
- Aggiornamento mappe statiche originali correggendo un solo settore per l albero non ancora tagliato
- Gestione generalizzata mancanza di SQI

* Wed Apr 12 2023 Daniele Branchini <dbranchini@arpae.it> - 0.19-1
- Fixed default fuzzypath

* Mon Apr  3 2023 Daniele Branchini <dbranchini@arpae.it> - 0.18-1
- Refactoring of `FUZZY_PATH` variable management
- Using TCLAP for input parameters in `stat_CleanID.cpp` and `RunCleanID.cpp`
- Fixing variable `check_undetect` (#23)

* Mon Mar 20 2023 Daniele Branchini <dbranchini@arpae.it> - 0.17-1
- Added data dir

* Mon Mar 13 2023 Daniele Branchini <dbranchini@arpae.it> - 0.16-1
- Update cleaning module with fuzzy logic methods
- Added `-U` flag in elaboradar to set bin as `undetect`
- Dropped CentOS7 support (#19)

* Wed Nov 30 2022 Daniele Branchini <dbranchini@arpae.it> - 0.15-1
- Added support for new proj API (#19)

* Fri Oct 30 2020 Daniele Branchini <dbranchini@arpae.it> - 0.14-1
- Fixed offset management
- Added AddCleanerQuantities
- Forcing clean to use undetect instead of nodata

* Thu Sep  3 2020 Daniele Branchini <dbranchini@arpae.it> - 0.13-1
- External reading for melting layer
- Flag Use_undetect

* Tue Apr 28 2020 Daniele Branchini <dbranchini@arpae.it> - 0.12-1
- Added RadarSite.h

* Wed Feb 19 2020 Daniele Branchini <dbranchini@arpae.it> - 0.11-1
- Removed melting layer top and bottom external reading top, now computed again

* Thu Mar 28 2019 Daniele Branchini <dbranchini@arpae.it> - 0.10-1
- Added undetect option to RunCleanID
- Added Ht and Hb parameters passage from rds

* Wed Sep 27 2017 Daniele Branchini <dbranchini@arpae.it> - 0.7-1
- Lots of stuff (please refer to upstream github logs)

* Wed Jul 13 2016 Daniele Branchini <dbranchini@arpae.it> - 0.6-1
- Bounding box now computed according to ODIM specification
- Closed #4
- Add RadarSite class to describe radar information site
- Some minor bug fixing

* Fri May 6 2016 Daniele Branchini <dbranchini@arpa.emr.it> - 0.5-1%{dist}
- Add method to save products with additional info in the filename
- Add test on presence of qual matrix before produce output product
- Add method to test cell_size
- Add method to write subimage with additional info in the filename
- Get RangeScale value from Odim scan

* Fri Apr 15 2016 Daniele Branchini <dbranchini@arpa.emr.it> - 0.4-1%{dist}
- Fixed bug in indexes calculation
- Fixed missing copy of h_radar attribute in resample_volume

* Tue Mar 22 2016 Daniele Branchini <dbranchini@arpa.emr.it> - 0.2-1%{dist}
- Closed #6, #7
- Fixed bug for range computation
- Removed dependency from site.h
- Changed samplign algorithm
- Double operation instead of unsigned operations
- Renamed lib
- Aadded doc

* Fri Jan 8 2016 Daniele Branchini <dbranchini@arpa.emr.it> - 0.1-1%{dist}
- First build
