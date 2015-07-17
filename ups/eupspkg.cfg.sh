build(){
    scons opt=3 prefix=$PREFIX version=$VERSION
}
install(){
    scons opt=3 install prefix=$PREFIX version=$VERSION
}
