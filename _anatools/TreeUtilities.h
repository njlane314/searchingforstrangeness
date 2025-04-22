#pragma once

#include "TTree.h"

namespace tree {
    template <typename T> class ManagedPointer : public std::unique_ptr<T> {
    public:
        ManagedPointer() : std::unique_ptr<T>(new T) {}

        T*& get_bare_ptr() {
            bare_ptr_ = this->get();
            return bare_ptr_;
        }

    protected:
        T* bare_ptr_ = nullptr;
    };

    template <typename T> void set_object_input_branch_address(TTree& in_tree, const std::string& branch_name, T*& address) {
        in_tree.SetBranchAddress(branch_name.c_str(), &address);
    }

    template <typename T> void set_object_input_branch_address(TTree& in_tree, const std::string& branch_name, ManagedPointer<T>& u_ptr) {
        T*& address = u_ptr.get_bare_ptr();
        set_object_input_branch_address(in_tree, branch_name, address);
    }

    void set_output_branch_address(TTree& out_tree, const std::string& branch_name, void* address, bool create = false, const std::string& leaf_spec = "") {
        if (create) {
            if (leaf_spec != "") {
                out_tree.Branch( branch_name.c_str(), address, leaf_spec.c_str() );
            } else {
                out_tree.Branch( branch_name.c_str(), address );
            }
        } else {
            out_tree.SetBranchAddress( branch_name.c_str(), address );
        }
    }

    template <typename T> void set_object_output_branch_address(TTree& out_tree, const std::string& branch_name, T*& address, bool create = false){
        if (create) 
            out_tree.Branch(branch_name.c_str(), &address);
        else 
            out_tree.SetBranchAddress(branch_name.c_str(), &address);
    }

    template <typename T> void set_object_output_branch_address(TTree& out_tree, const std::string& branch_name, ManagedPointer<T>& u_ptr, bool create = false) {
        T*& address = u_ptr.get_bare_ptr();
        set_object_output_branch_address(out_tree, branch_name, address, create);
    }
}