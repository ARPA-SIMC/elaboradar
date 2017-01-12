Summary:	Library and tools to handle weather radar images and data
Name: 		elaboradar
Version: 	0.7
Release: 	1
License: 	GPL
Group: 		Applications/Meteo
URL:            https://github.com/arpa-simc/%{name}
Source0:        https://github.com/arpa-simc/%{name}/archive/v%{version}-%{release}.tar.gz#/%{name}-%{version}-%{release}.tar.gz
BuildRoot:      %{_tmppath}/%{name}-%{version}-%{release}-buildroot
Packager: 	Daniele Branchini <dbranchini@arpa.emr.it>
BuildRequires:	gcc-c++, hdf5-devel, eigen3-devel, radarlib-devel, log4c-devel, gsl-devel, tclap-devel
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
%setup -q -n %{name}-%{version}-%{release}
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
%{_bindir}/classificatore
%{_bindir}/elaboradar

%changelog
* Tue Jan 12 2017 Daniele Branchini <dbranchini@arpae.it> - 0.7-1
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
