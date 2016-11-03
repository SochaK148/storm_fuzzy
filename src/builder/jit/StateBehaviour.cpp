#include "src/builder/jit/StateBehaviour.h"

#include "src/adapters/CarlAdapter.h"

namespace storm {
    namespace builder {
        namespace jit {
            
            template <typename IndexType, typename ValueType>
            StateBehaviour<IndexType, ValueType>::StateBehaviour() : expanded(false) {
                // Intentionally left empty.
            }
            
            template <typename IndexType, typename ValueType>
            void StateBehaviour<IndexType, ValueType>::addChoice(Choice<IndexType, ValueType>&& choice) {
                choices.emplace_back(std::move(choice));
            }
            
            template <typename IndexType, typename ValueType>
            Choice<IndexType, ValueType>& StateBehaviour<IndexType, ValueType>::addChoice() {
                choices.emplace_back();
                return choices.back();
            }
            
            template <typename IndexType, typename ValueType>
            typename StateBehaviour<IndexType, ValueType>::ContainerType const& StateBehaviour<IndexType, ValueType>::getChoices() const {
                return choices;
            }
            
            template <typename IndexType, typename ValueType>
            void StateBehaviour<IndexType, ValueType>::addStateReward(ValueType const& stateReward) {
                stateRewards.push_back(stateReward);
            }
            
            template <typename IndexType, typename ValueType>
            void StateBehaviour<IndexType, ValueType>::addStateRewards(std::vector<ValueType>&& stateRewards) {
                this->stateRewards = std::move(stateRewards);
            }

            template <typename IndexType, typename ValueType>
            std::vector<ValueType> const& StateBehaviour<IndexType, ValueType>::getStateRewards() const {
                return stateRewards;
            }
            
            template <typename IndexType, typename ValueType>
            void StateBehaviour<IndexType, ValueType>::reduce(storm::jani::ModelType const& modelType) {
                if (choices.size() > 1) {
                    if (modelType == storm::jani::ModelType::DTMC || modelType == storm::jani::ModelType::CTMC) {
                        std::size_t totalCount = choices.size();
                        for (auto it = ++choices.begin(), ite = choices.end(); it != ite; ++it) {
                            choices.front().add(std::move(*it));
                        }
                        choices.resize(1);
                        choices.front().compress();

                        if (modelType == storm::jani::ModelType::DTMC) {
                            choices.front().divideDistribution(static_cast<ValueType>(totalCount));
                        }
                    } else {
                        for (auto& choice : choices) {
                            choice.compress();
                        }
                    }
                } else if (choices.size() == 1) {
                    choices.front().compress();
                }
            }
            
            template <typename IndexType, typename ValueType>
            bool StateBehaviour<IndexType, ValueType>::isExpanded() const {
                return expanded;
            }
            
            template <typename IndexType, typename ValueType>
            void StateBehaviour<IndexType, ValueType>::setExpanded() {
                expanded = true;
            }
            
            template <typename IndexType, typename ValueType>
            bool StateBehaviour<IndexType, ValueType>::empty() const {
                return choices.empty();
            }
            
            template <typename IndexType, typename ValueType>
            std::size_t StateBehaviour<IndexType, ValueType>::size() const {
                return choices.size();
            }

            template <typename IndexType, typename ValueType>
            void StateBehaviour<IndexType, ValueType>::clear() {
                choices.clear();
                expanded = false;
            }
            
            template class StateBehaviour<uint32_t, double>;
            template class StateBehaviour<uint32_t, storm::RationalNumber>;
            template class StateBehaviour<uint32_t, storm::RationalFunction>;
            
        }
    }
}
