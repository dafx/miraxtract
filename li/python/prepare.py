from scipy import signal



def array_to_c(arr, name):
    print(f"\nconst float {name}[] = {{")
    for i in range(len(arr) - 1):
        print(f"    /*{i:02d}*/ {arr[i]:e}f,")
    print(f"    /*{len(arr) - 1:02d}*/ {arr[-1]:e}f")    
    print("};\n")


def gen_filter_coefs():
    sos = signal.butter(8, 0.2, output='sos')
    array_to_c(sos.flatten(), "deci_4_filter_coefs")
    zi = signal.sosfilt_zi(sos)
    array_to_c(zi.flatten(), "deci_4_filter_zi")
    return sos


if __name__ == '__main__':
    gen_filter_coefs()
