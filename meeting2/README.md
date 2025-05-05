# Αναφορά profiling με Vtune

## 1. Top-level functions
Αρχικά αναφέρονται οι βασικές και πιο χρονοβόρες top-level συναρτήσεις των applications που έτρεξα. Είναι λογικό να φαίνονται στα reports ως πολύ χρονοβόρες, καθώς καλούν μια σειρά από bottom level συναρτήσεις για να κάνουν τις πράξεις. Στην ks_imp_rhmc εφαρμογή, είτε στα HISQ builds (πιο αργά αλλά πιο accurate) είτε στα πιο απλά builds, οι top level συναρτήσεις που δεν είναι κοινές, καλούσαν τις ίδιες bottom level συναρτήσεις ούτως ή άλλως. Στην ks_imp_dyn εφαρμογή, οι συναρτήσεις top level ήταν είτε οι ίδιες με αυτές της ks_imp_rhmc, είτε ελάχιστα διαφορετικά versions αυτών, που πάλι κατά βάση καλούσαν τις ίδιες bottom level συναρτήσεις. 
Έχω τον παρακάτω πίνακα:

| Function | Τι κάνει | Τι bottom level καλεί |
|----------|-------|-------|
| dslash_fn_field_special | Εφαρμόζει τον Dirac operator σε πεδίο fermion. | mult_su3_mat_vec_sum_4dir, mult_su3_mat_vec, add_su3_vector, sub_su3_vector |
| fermion_force_fn_multi | Υπολογίζει τη συμβολή των fermions στην gauge δύναμη. | su3_projector, mult_su3_na, scalar_mult_add_su3_matrix|
| path_product_fields | Υπολογίζει το γινόμενο SU(3) συνδέσμων κατά μήκος μονοπατιού. | mult_su3_nn, mult_su3_na |
| compute_gen_staple_field | Κατασκευάζει το staple (3-link term) για βελτιωμένες δράσεις. | mult_su3_na, mult_su3_nn, scalar_mult_add_su3_matrix |
| fn_fermion_force_multi_hisq_smearing | HISQ έκδοση του fermion force με smearing και Naik όρους. | mult_su3_na, scalar_mult_add_su3_matrix, su3_projector |
| link_transport_connection | Μεταφέρει SU(3) πίνακα κατά μία κατεύθυνση. | mult_su3_nn, mult_su3_an |
| link_gather_connection_hisq | Μαζεύει SU(3) από γειτονικό site. | Δεν καλεί κάποια computationally heavy bottom level συνάρτηση, παρόλο που κάποιες φορές φαίνεται να παίρνει χρόνο |

---

## 2. Bottom-level functions
Εδώ φαίνονται οι πιο βασικές bottom level συναρτήσεις που εντώπισα, είτε αυτές φαίνονται άμεσα στα vtune reports, είτε εμφανίζονται έμμεσα στις top level. Την υπολογιστική βαρύτητα τη σκέφτηκα με βάση τα reports, έναν πρόχειρο υπολογισμό από FLOPS σε κάθε συνάρτηση και κάπως διαισθητικά, επομένως μπορεί να μην είναι 100% ορθή.

| Συνάρτηση | Περιγραφή | Υπολογιστική Βαρύτητα |
|-----------|-----------|----------|
| [mult_su3_nn](https://github.com/JiMan5/HW_SW_codesign_thesis/blob/main/bottom_lvl_code/m_mat_nn.c) | Κανονικό γινόμενο SU(3) πινάκων. | Πολύ βαριά |
| [mult_su3_na](https://github.com/JiMan5/HW_SW_codesign_thesis/blob/main/bottom_lvl_code/m_mat_na.c) | Γινόμενο SU(3) πινάκων, πρώτος κανονικός, δεύτερος adjoint | Πολύ βαριά |
| [mult_su3_an](https://github.com/JiMan5/HW_SW_codesign_thesis/blob/main/bottom_lvl_code/m_mat_an.c) | Γινόμενο SU(3) πινάκων, πρώτος adjoint, δεύτερος κανονικός | Πολύ βαριά |
| [mult_su3_mat_vec_sum_4dir](https://github.com/JiMan5/HW_SW_codesign_thesis/blob/main/bottom_lvl_code/m_mv_s_4dir.c) | Πολλαπλασιάζει 4 SU(3) πίνακες με διανύσματα και προσθέτει τα παραγώμενα διανύσματα σε ένα τελικό διάνυσμα | Βαριά |
| [mult_su3_mat_vec](https://github.com/JiMan5/HW_SW_codesign_thesis/blob/main/bottom_lvl_code/m_matvec.c) | Πολλαπλασιασμός πίνακα SU(3) με διάνυσμα su3_vector | Μέτρια | 
| [mult_adj_su3_mat_vec](https://github.com/JiMan5/HW_SW_codesign_thesis/blob/main/bottom_lvl_code/m_amatvec.c) | Πολλαπλασιασμός adjoint πίνακα SU(3) με διάνυσμα su3_vector | Μέτρια |
| [su3_projector](https://github.com/JiMan5/HW_SW_codesign_thesis/blob/main/bottom_lvl_code/su3_proj.c) | Εξωτερικό γινόμενο διανυσμάτων su3_vector | Μέτρια |
| [scalar_mult_add_su3_matrix](https://github.com/JiMan5/HW_SW_codesign_thesis/blob/main/bottom_lvl_code/s_m_a_mat.c) | Πολλαπλασιασμός αριθμού με πίνακα και έπειτα πρόσθεση πινάκων. (C = A + s*B, A,B --> SU(3)| Μέτρια |
| [scalar_mult_sub_su3_vector](https://github.com/JiMan5/HW_SW_codesign_thesis/blob/main/bottom_lvl_code/s_m_s_vec.c) | Πολλαπλασιασμός αριθμού με διάνυσμα και έπειτα πρόσθεση διανυσμάτων. (C = A + s*B, A,B --> su3_vector| Μέτρια |
| [add_su3_matrix](https://github.com/JiMan5/HW_SW_codesign_thesis/blob/main/bottom_lvl_code/addmat.c) | Άθροιση πινάκων. | Ελαφριά |
| [sub_su3_matrix](https://github.com/JiMan5/HW_SW_codesign_thesis/blob/main/bottom_lvl_code/submat.c) | Αφαίρεση πινάκων. | Ελαφριά |
| [add_su3_vector](https://github.com/JiMan5/HW_SW_codesign_thesis/blob/main/bottom_lvl_code/addvec.c) | Άθροιση διανυσμάτων. | Ελαφριά |
| [sub_su3_vector](https://github.com/JiMan5/HW_SW_codesign_thesis/blob/main/bottom_lvl_code/subvec.c) | Αφαίρεση διανυσμάτων. | Ελαφριά | 
| [su3_adjoint](https://github.com/JiMan5/HW_SW_codesign_thesis/blob/main/bottom_lvl_code/su3_adjoint.c) | Adjoint SU(3) πίνακα. | Ελαφριά |
| [su3mat_copy](https://github.com/JiMan5/HW_SW_codesign_thesis/blob/main/bottom_lvl_code/su3mat_copy.c) | Αντιγραφή SU(3) πίνακα. | Ελαφριά |

---

## 3. Half-Wilson bottom-level functions

Τέλος, στην εφαρμογή ks_imp_dyn, στα reports φαίνονται συναρτήσεις για πράξεις μεταξύ SU(3) πινάκων και διανυσμάτων, οι οποίες όμως χρησιμοποιούν τα optimized Half-Wilson vectors. Αυτές οι συναρτήσεις παίρνουν πολύ χρόνο, όπως φαίνεται στα αντίστοιχα reports και είναι όλες bottom level και computationally heavy, επομένως θεωρητικά θα μπορούσαν να επιταχυνθούν σε FPGA kernels. Όμως, στην πράξη, λόγω του optimization που ήδη χρησιμοποιόυν, δεν είμαι σίγουρος για το αν είναι καλοί υποψήφιοι για HW kernels, θα ήθελα τη γνώμη σας πάνω σε αυτές. Η πρώτη μου σκέψη είναι πως δε θα είχαμε κάποιο θέμα αφού αποτελούν ουσιαστικά hand-unrolled loops, πράγμα που σημαίνει πως μπορούμε να τις παραλληλοποιήσουμε. Όλες οι συναρτήσεις αυτού του τύπου βρίσκονται στο αρχείο [ff_opt.c](https://github.com/JiMan5/HW_SW_codesign_thesis/blob/main/bottom_lvl_code/ff_opt.c), όμως μόνο κάποιες από αυτές εμφανίζονται άμεσα στα reports.

## 4. Σκέψεις προς συζήτηση/έρευνα 

Φυσικά, παρατηρώ πως οι bottom level συναρτήσεις είναι εξ αρχής μικρές και πολύ εύκολα παραλληλοποιησιμες. Είτε μιλάμε για SU(3)*SU(3),ειτε μιλάμε για SU(3)*vector, είναι λίγα τα FLOPS σε κάθε κάλεσμα των συναρτήσεων αυτών. Παρόλα αυτά, αυτές οι συναρτήσεις καλούνται πολλές χιλιάδες/εκατομμύρια φορές από τις top level functions, ανάλογα με το μέγεθος του lattice. Για αυτό, αξίζει να γίνει μια έρευνα πίσω από τα data dependencies κλπ, ώστε loops που καλούν αυτες τις bottom level συναρτησεις να γίνουν pipelined επισης, παρόλο που βρίσκονται εντός των top level συναρτήσεων. Φαντάζομαι εντός του SoC που θα υλοποιηθεί στα πλαίσια της διπλωματικής εργασίας, θα μελετηθουν τέτοια κομμάτια, data flow, dependencies κα. 
