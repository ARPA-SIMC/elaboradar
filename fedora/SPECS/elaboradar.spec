Summary:	Library and tools to handle weather radar images and data
Name: 		elaboradar
Version: 	0.2
Release: 	1
License: 	GPL
Group: 		Applications/Meteo
URL:            https://github.com/arpa-simc/%{name}
Source0:        https://github.com/arpa-simc/%{name}/archive/v%{version}-%{release}.tar.gz#/%{name}-%{version}-%{release}.tar.gz
BuildRoot:      %{_tmppath}/%{name}-%{version}-%{release}-buildroot
Packager: 	Daniele Branchini <dbranchini@arpa.emr.it>
BuildRequires:	gcc-c++, hdf5-devel, eigen3-devel, radarlib-devel
Requires:       hdf5

%description
Library and to handle weather radar images and data

%package devel
Requires: elaboradar = %{version}
Group: Libraries/Meteo
Summary: Development for radarelab library

%description devel
Development for elaboradar library

#package doc
#Summary: elaboradar documentation
#Group: Libraries/Meteo

#description doc
#elaboradar library documentation

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
%{_libdir}/lib%{name}.so.0*

%files devel
%defattr(-,root,root,-)
%{_includedir}/%{name}/*
%{_libdir}/lib%{name}.a
%{_libdir}/lib%{name}.la
%{_libdir}/lib%{name}.so
%{_libdir}/pkgconfig/%{name}.pc

#files doc
#defattr(-,root,root,-)
#doc %{_docdir}/%{name}

%files tools
%{_bindir}/classificatore
%{_bindir}/elaboradar

%changelog
* Tue Mar 22 2016 Daniele Branchini <dbranchini@arpa.emr.it> - 0.2-1%{dist}
- closed \#6, \#7
- fixed bug for range computation
- removed dependency from site.h
- changed samplign algorithm
- double operation instead of unsigned operations

* Fri Jan 8 2016 Daniele Branchini <dbranchini@arpa.emr.it> - 0.1-1%{dist}
- First build
