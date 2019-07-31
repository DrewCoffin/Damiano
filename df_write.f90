        SUBROUTINE df_write(index)
        use indeces
        use df
        implicit none
        character(len=7) :: fname
        integer :: index
        call getname('d_df',index,fname)
        open(10,file=fname,form='unformatted',status='unknown')
        write(10) df_array
        close(10)
        return
        END SUBROUTINE df_write

