# Note: define srcarchivename in CI build only.
%{!?srcarchivename: %global srcarchivename %{name}-%{version}-%{release}}

Summary:	Library and tools to handle weather radar images and data
Name: 		elaboradar
Version: 	0.14
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
%{_libdir}/libradarelab.la
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

%changelog
* Fri Oct 30 2020 Daniele Branchini <dbranchini@arpae.it> - 0.14-1
- fixed offset management
- added AddCleanerQuantities
- forcing clean to use undetect instead of nodata

* Thu Sep  3 2020 Daniele Branchini <dbranchini@arpae.it> - 0.13-1
- external reading for melting layer
- flag Use_undetect

* Tue Apr 28 2020 Daniele Branchini <dbranchini@arpae.it> - 0.12-1
- added RadarSite.h

* Wed Feb 19 2020 Daniele Branchini <dbranchini@arpae.it> - 0.11-1
- removed melting layer top and bottom external reading top, now computed again

* Thu Mar 28 2019 Daniele Branchini <dbranchini@arpae.it> - 0.10-1
- added undetect option to RunCleanID
- added Ht and Hb parameters passage from rds

* Wed Sep 27 2017 Daniele Branchini <dbranchini@arpae.it> - 0.7-1
- lots of stuff (please refer to upstream github logs)

* Wed Jul 13 2016 Daniele Branchini <dbranchini@arpae.it> - 0.6-1
- bounding box now computed according to ODIM specification
- closed #4
- add RadarSite class to describe radar information site
- some minor bug fixing

* Fri May 6 2016 Daniele Branchini <dbranchini@arpa.emr.it> - 0.5-1%{dist}
- Add method to save products with additional info in the filename
- Add test on presence of qual matrix before produce output product
- Add method to test cell_size
- Add method to write subimage with additional info in the filename
- Get RangeScale value from Odim scan

* Fri Apr 15 2016 Daniele Branchini <dbranchini@arpa.emr.it> - 0.4-1%{dist}
- fixed bug in indexes calculation
- fixed missing copy of h_radar attribute in resample_volume

* Tue Mar 22 2016 Daniele Branchini <dbranchini@arpa.emr.it> - 0.2-1%{dist}
- closed \#6, \#7
- fixed bug for range computation
- removed dependency from site.h
- changed samplign algorithm
- double operation instead of unsigned operations
- renamed lib
- added doc

* Fri Jan 8 2016 Daniele Branchini <dbranchini@arpa.emr.it> - 0.1-1%{dist}
- First build
